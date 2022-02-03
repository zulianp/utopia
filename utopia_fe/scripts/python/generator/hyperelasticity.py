# import sympy as sp
from sympy import *
from sympy.utilities.codegen import codegen
from sympy.utilities.codegen import CCodeGen
from sympy.codegen.ast import Assignment
from sympy.codegen.ast import AddAugmentedAssignment
from sympy.printing.c import C99CodePrinter
# from sympy.codegen.cfunctions import *
import os
import sys

from utopia_mech import *


class HyperElasticModel:
	def __init__(self, d):
		self.d = d
		self.F = deformation_gradient(d)
		self.J = det(self.F)
		self.F_inv = Inverse(self.F)
		self.F_inv_t = self.F_inv.T
		self.C = self.F.T*self.F
		self.I_C = trace(self.C)

		self.strain_energy_function = 0

###############################
# NeoHookean models
###############################

class NeoHookeanOgden(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		I_C = self.I_C
		J = self.J

		mu, lmbda = symbols('mu lambda')
		self.params = [(mu, 1.0), (lmbda,1.0)]

		self.fun = mu/2 *(I_C - d) - mu * log(J) + (lmbda/2) * (log(J))**2
		self.name = 'NeoHookeanOgden'


class NeoHookeanBower(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		I_C = self.I_C
		J = self.J

		mu, lmbda = symbols('mu lambda')
		self.params = [(mu, 1.0), (lmbda, 1.0)]

		self.fun = mu/2 *(J**(-Rational(2, 3))*I_C - d) + (lmbda/2) * (J-1)**2
		self.name = 'NeoHookeanBower'


class NeoHookeanWang(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		I_C = self.I_C
		J = self.J

		mu, lmbda = symbols('mu lambda')
		self.params = [(mu, 1.0), (lmbda, 1.0)]

		self.fun = mu/2 *(J**(-Rational(2, 3))*I_C - d) + (lmbda/2) * (J-1) 
		self.name = 'NeoHookeanWang'


class NeoHookeanSmith(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		I_C = self.I_C
		J = self.J

		mu, lmbda = symbols('mu lambda')
		self.params = [(mu, 1.0), (lmbda, 1.0)]

		mu_s = Rational(4,3)*mu
		lmbda_s = lmbda + Rational(5,6)*lmbda
		alpha = 1 + mu_s/lmbda_s - mu_s/(4*lmbda_s)

		self.fun = mu_s/2 * (I_C - d) + lmbda_s/2 * (J - alpha)**2 - mu_s/2 * log(I_C + 1)
		self.name = 'NeoHookeanSmith'

###############################
# Fung (https://en.wikipedia.org/wiki/Soft_tissue#Fung-elastic_material)
###############################
class Fung(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		I_C = self.I_C
		J = self.J

		a, b, c = symbols('a b c')
		self.params = [(a, 1.0), (b, 1.0), (c, 1.0)]

		self.fun = Rational(1,2) * (a * (I_C - d) + b * (exp(c*(I_C - d) - 1)))
		self.name = 'Fung'


###############################
# Mooney-Rivlin (https://en.wikipedia.org/wiki/Mooney%E2%80%93Rivlin_solid)
###############################
class MooneyRivlin(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		I_C = self.I_C
		J = self.J
		C = self.C

		C1, C2 = symbols('C1 C2')
		self.params = [(C1, 1.0), (C2, 1.0)]
		
		CikxCki = 0
		for i in range(0, d):
			for k in range(0, d):
				CikxCki += C[i,k] * C[k,i]

		I_C2 = Rational(1, 2) * (I_C**2 - CikxCki)
		# I_C3 = det(C)

		I1 = J**(-Rational(d-1, d))*I_C
		I2 = J**(-Rational(d+1, d))*I_C2;
		self.fun = C1 * (I1 - d) + C2 * (I2 - d)
		self.name = 'MooneyRivlin'


class HyperElasticity:

	# def __init__(self):

	def generate_files(self, model, output_dir, simplify_expressions):
		# unpack
		params = model.params
		W = model.fun 
		F = model.fun 
		name = model.name
		F = model.F

		#############################################
		# UI
		#############################################

		if not os.path.exists(output_dir):
			print(f"Creating directory {output_dir}")
			os.mkdir(output_dir)


		#############################################
		# FE
		#############################################

		d = F.shape[0]
		grad_trial = trial_gradient(d)
		grad_test = test_gradient(d)

		dx = symbols('dx')
		
		
		#############################################
		# Energy
		#############################################
		
		energy = W * dx

		#############################################
		# First-Piola-Kirchoff stress tensor
		#############################################
		
		P = first_piola(W, F)

		#############################################
		# Gradient
		#############################################
		print("Creating Gradient")

		linear_form = 0
		for i in range(0, d):
			for j in range(0, d):
		 		linear_form += P[i,j] * grad_test[i,j]

		linear_form *= dx

		#############################################
		# Hessian
		#############################################
		print("Creating Hessian")

		contraction = 0
		for i in range(0, d):
			for j in range(0, d):
		 		contraction += P[i,j] * grad_trial[i,j]

		bilinear_form = 0

		for i in range(0, d):
			for j in range(0, d):
				Hij = diff(contraction, F[i, j]);
				bilinear_form += Hij * grad_test[i,j]

		bilinear_form *= dx
		#############################################
		# Substitutions
		#############################################
		print("Substituting variables")

		tp = TensorProductBasis(d)

		hessian_expression_list = []
		gradient_expression_list = []
		energy_expression_list = []



		for d1 in range(0, d):
			subsituted = tp.linear_subs("test", d1, linear_form)

			if simplify_expressions:
				subsituted = simplify(subsituted)

			gradient_expression_list.append(AddAugmentedAssignment(symbols(f"lf[{d1}]"), subsituted))

			for d2 in range(0, d):
				subsituted = tp.bilinear_subs("test", d1, "trial", d2, bilinear_form)

				if simplify_expressions:
					subsituted = simplify(subsituted)

				hessian_expression_list.append(AddAugmentedAssignment(symbols(f"bf[{d1*d + d2}]"), subsituted))


		energy_expression_list.append(AddAugmentedAssignment(symbols("e"), energy))


		full_expression_list = []
		full_expression_list.extend(hessian_expression_list)
		full_expression_list.extend(gradient_expression_list)
		full_expression_list.extend(energy_expression_list)


		#############################################
		# Generate code
		#############################################

		print("Generating code")

		

		generator = KernelGenerator(d)
		generator.generate_class(
			Template(
				"templates/utopia_tpl_hyperelasticity.hpp",
				"templates/utopia_tpl_hyperelasticity_impl.hpp",
				f"{output_dir}/../utopia_hyperelasticity_{name}.hpp",
				f"{output_dir}/utopia_hyperelasticity_{name}_{d}.hpp"), 
			name, 
			params,
	        energy_expression_list,
	        gradient_expression_list,
	        hessian_expression_list,
	        output_dir)

		# generator.generate_files(hessian_expression_list, 'templates/utopia_tpl_elasticity_hessian.hpp', output_dir + '/generated_AutoHyperElasticityHessian.hpp')
		# generator.generate_files(gradient_expression_list, 'templates/utopia_tpl_elasticity_gradient.hpp', output_dir + '/generated_AutoHyperElasticityGradient.hpp')
		# generator.generate_files(energy_expression_list, 'templates/utopia_tpl_elasticity_energy.hpp', output_dir + '/generated_AutoHyperElasticityEnergy.hpp')
		# generator.generate_files(full_expression_list, 'templates/utopia_tpl_elasticity.hpp', output_dir + '/generated_AutoHyperElasticity.hpp')


def generate_materials(d,simplify_expressions):
	output_dir = f'../../../backend/kokkos/assembly/mech/generated/{d}D'
	# models = [NeoHookeanOgden(d), NeoHookeanBower(d), NeoHookeanWang(d), NeoHookeanSmith(d), Fung(d), MooneyRivlin(d)]
	models = [NeoHookeanOgden(d)]

	for m in models:
		HyperElasticity().generate_files(m, output_dir, simplify_expressions)	

def main(args):
	generate_materials(2)
	generate_materials(3)

if __name__ == '__main__':
    main(sys.argv[1:])
