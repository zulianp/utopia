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
		self.block_size = d
		self.F = deformation_gradient(d)
		self.J = det(self.F)
		self.F_inv = Inverse(self.F)
		self.F_inv_t = self.F_inv.T
		self.C = self.F.T*self.F
		self.I_C = trace(self.C)
		self.fun = 0

		self.form = zeros(3)
		self.use_default_parameter_reader = False

	def compute_forms(self):
		params = self.params
		W = self.fun
		name = self.name
		F = self.F
		d = self.d

		#############################################
		# FE
		#############################################

		d = F.shape[0]
		grad_trial = trial_gradient(d)
		grad_test = test_gradient(d)

		dx = symbols('dx')
		self.dx = dx


		#############################################
		# Energy
		#############################################

		energy = W * dx

		self.form[0] = energy
		#############################################
		# First-Piola-Kirchoff stress tensor
		#############################################

		P = first_piola(W, F)
		self.P = P

		#############################################
		# Gradient
		#############################################
		print("Creating Gradient")

		linear_form = 0
		for i in range(0, d):
			for j in range(0, d):
		 		linear_form += P[i,j] * grad_test[i,j]

		linear_form *= dx

		self.form[1] = linear_form

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

		self.form[2] = bilinear_form

class IncompressibleHyperElasticModel(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)
		self.block_size = d + 1 # Add one element for the pressure
		self.p = symbols('p')
		# self.J = symbols('J')

	def compute_forms(self):
		super().compute_forms()

		params = self.params
		W = self.fun
		name = self.name
		J = self.J
		d = self.d
		p = self.p
		dx = self.dx
		trial = trial_function(1)
		test = test_function(1)

		#############################################
		# Energy
		#############################################

		# should be complete with subclass model

		print(self.fun)

		#############################################
		# Gradient
		#############################################
		print("Gradient(pressure)")

		# dWdJ = simplify(diff(self.fun, J))
		dWdp = simplify(diff(self.fun, p)) * test * dx
		# print(self.form[1])
		print(dWdp)

		# self.form[1] += (dWdJ)

		#############################################
		# Hessian
		#############################################
		print("Hessian(pressure)")
		dWdp2 = diff(dWdp, p)
		print(dWdp2)
		# self.form[1] += (dWdJ)


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
		self.use_default_parameter_reader = True


class NeoHookeanBower(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		I_C = self.I_C
		J = self.J

		mu, lmbda = symbols('mu lambda')
		self.params = [(mu, 1.0), (lmbda, 1.0)]

		self.fun = mu/2 *(J**(-Rational(2, 3))*I_C - d) + (lmbda/2) * (J-1)**2
		self.name = 'NeoHookeanBower'
		self.use_default_parameter_reader = True


class NeoHookeanWang(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		I_C = self.I_C
		J = self.J

		mu, lmbda = symbols('mu lambda')
		self.params = [(mu, 1.0), (lmbda, 1.0)]

		self.fun = mu/2 *(J**(-Rational(2, 3))*I_C - d) + (lmbda/2) * (J-1) 
		self.name = 'NeoHookeanWang'
		self.use_default_parameter_reader = True


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
		self.use_default_parameter_reader = True

###############################
# Fung (https://en.wikipedia.org/wiki/Soft_tissue#Fung-elastic_material)
###############################
class Fung(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		a, b, c, k = symbols('a b c k')
		self.params = [(a, 1.0), (b, 1.0), (c, 1.0), (k, 1)]

		C = self.C
		J = self.J
		E = (C - eye(d, d))/2

		# Isotropic
		trE = trace(E)
		q = a * trE
		Q = b * trE

		self.fun = Rational(1,2) * (q + c * (exp(Q) - 1)) + (k/2)*(J - 1)**2
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


class IncompressibleMooneyRivlin(IncompressibleHyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		I_C = self.I_C
		J = self.J
		C = self.C
		p = self.p

		C1, C2 = symbols('C1 C2')
		self.params = [(C1, 1.0), (C2, 1.0)]

		CikxCki = 0
		for i in range(0, d):
			for k in range(0, d):
				CikxCki += C[i,k] * C[k,i]

		I_C2 = Rational(1, 2) * (I_C**2 - CikxCki)

		self.fun = C1 * (I_C - d) + C2 * (I_C2 - d) + p/2 * (J - 1)**2
		self.name = 'IncompressibleMooneyRivlin'

class HyperElasticity:

	# def __init__(self):

	def generate_files(self, model, output_dir, simplify_expressions):
		#############################################
		# UI
		#############################################

		if not os.path.exists(output_dir):
			print(f"Creating directory {output_dir}")
			os.mkdir(output_dir)

		#############################################
		# FE
		#############################################

		model.compute_forms()
		tp = TensorProductBasis(model.block_size)

		hessian_expression_list = []
		gradient_expression_list = []
		energy_expression_list = []

		for d1 in range(0, model.block_size):
			subsituted = tp.linear_subs("test", d1, model.form[1])

			if simplify_expressions:
				subsituted = simplify(subsituted)

			gradient_expression_list.append(AddAugmentedAssignment(symbols(f"lf[{d1}]"), subsituted))

			for d2 in range(0, model.block_size):
				subsituted = tp.bilinear_subs("test", d1, "trial", d2, model.form[2])

				if simplify_expressions:
					subsituted = simplify(subsituted)

				hessian_expression_list.append(AddAugmentedAssignment(symbols(f"bf[{d1*model.block_size + d2}]"), subsituted))


		energy_expression_list.append(AddAugmentedAssignment(symbols("e"), model.form[0]))

		full_expression_list = []
		full_expression_list.extend(hessian_expression_list)
		full_expression_list.extend(gradient_expression_list)
		full_expression_list.extend(energy_expression_list)

		#############################################
		# Generate code
		#############################################

		print("Generating code")

		generator = KernelGenerator(model.d)
		generator.generate_class(
			Template(
				"templates/utopia_tpl_hyperelasticity.hpp",
				"templates/utopia_tpl_hyperelasticity_impl.hpp",
				f"{output_dir}/../utopia_hyperelasticity_{model.name}.hpp",
				f"{output_dir}/utopia_hyperelasticity_{model.name}_{model.d}.hpp"),
			model,
	        energy_expression_list,
	        gradient_expression_list,
	        hessian_expression_list,
	        output_dir)

def generate_materials(d,simplify_expressions):
	output_dir = f'../../../backend/kokkos/assembly/mech/generated/{d}D'
	# models = [NeoHookeanOgden(d), NeoHookeanBower(d), NeoHookeanWang(d), NeoHookeanSmith(d), Fung(d), MooneyRivlin(d)]
	# models = [NeoHookeanOgden(d)]
	# models = [Fung(d)]
	models=[IncompressibleMooneyRivlin(d)]

	for m in models:
		HyperElasticity().generate_files(m, output_dir, simplify_expressions)	

def main(args):
	generate_materials(2,True)
	# generate_materials(3, False)

if __name__ == '__main__':
    main(sys.argv[1:])
