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


# class HyperElasticModel:
# 	def __init__(self, d):
# 		self.d = d
# 		self.F = deformation_gradient(d)
# 		self.J = det(self.F)
# 		self.F_inv = Inverse(self.F)
# 		self.F_inv_t = self.F_inv.T
# 		self.C = self.F.T*self.F
# 		self.I_C = trace(self.C)

# 		self.strain_energy_function = 0


class HyperElasticity:

	# def __init__(self):

	def generate_files(self, params, W, F, name, output_dir, simplify_expressions):
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
				# Hij = simplify(diff(contraction, F[i, j]));
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



def main(args):
	d = 3

	if d==2:
		output_dir = '../../../backend/kokkos/assembly/mech/generated/2D'
	else:
		output_dir = '../../../backend/kokkos/assembly/mech/generated/3D'


	###############################
	# Nonlinear quantities
	###############################
	F = deformation_gradient(d)
	J = det(F)
	F_inv = Inverse(F)
	F_inv_t = F_inv.T
	C = F.T*F
	I_C = trace(C)

	###############################
	# Material parameters
	###############################

	mu, lmbda = symbols('mu lambda')
	
	###############################
	# Strain energy function
	###############################

	###############################
	# NeoHookean models
	###############################

	Ogden = mu/2 *(I_C - d) - mu * log(J) + (lmbda/2) * (log(J))**2
	Bower = mu/2 *(J**(-Rational(2, 3))*I_C - d) + (lmbda/2) * (J-1)**2
	Wang  = mu/2 *(J**(-Rational(2, 3))*I_C - d) + (lmbda/2) * (J-1) # Does not work

	mu_s = Rational(4,3)*mu
	lmbda_s = lmbda + Rational(5,6)*lmbda
	alpha = 1 + mu_s/lmbda_s - mu_s/(4*lmbda_s)

	Smith = mu_s/2 * (I_C - d) + lmbda_s/2 * (J - alpha)**2 - mu_s/2 * log(I_C + 1)

	###############################
	# Fung (https://en.wikipedia.org/wiki/Soft_tissue#Fung-elastic_material)
	###############################

	a, b, c = symbols('a b c')
	Fung = Rational(1,2) * (a * (I_C - d) + b * (exp(c*(I_C - d) - 1)))

	###############################
	# Mooney-Rivlin (https://en.wikipedia.org/wiki/Mooney%E2%80%93Rivlin_solid)
	###############################

	C1, C2 = symbols('C1 C2')
	
	CikxCki = 0
	for i in range(0, d):
		for k in range(0, d):
			CikxCki += C[i,k] * C[k,i]

	I_C2 = Rational(1, 2) * (I_C**2 - CikxCki)
	# I_C3 = det(C)

	I1 = J**(-Rational(d-1, d))*I_C
	I2 = J**(-Rational(d+1, d))*I_C2;
	MooneyRivlin = C1 * (I1 - d) + C2 * (I2 - d)

	###############################
	# Select model
	###############################

	strain_energy_function = MooneyRivlin
	name = "MooneyRivlin"
	params = [C1, C2]

	# strain_energy_function = Fung
	# name = "Fung"
	# params = [a,b,c]

	# strain_energy_function = Ogden
	# name = "NeohookeanOgden"
	# params = [mu, lmbda]

	# strain_energy_function = Bower
	# name = "NeohookeanBower"
	# params = [mu, lmbda]

	# strain_energy_function = Wang
	# name = "NeoHookeanWang"
	# params = [mu, lmbda]

	# strain_energy_function = Smith
	# name = "NeohookeanSmith"
	# params = [mu, lmbda]

	HyperElasticity().generate_files(params, strain_energy_function, F, name, output_dir, False)

if __name__ == '__main__':
    main(sys.argv[1:])
