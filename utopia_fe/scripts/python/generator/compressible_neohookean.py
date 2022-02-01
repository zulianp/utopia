# import sympy as sp
from sympy import *
from sympy.utilities.codegen import codegen
from sympy.utilities.codegen import CCodeGen
from sympy.codegen.ast import Assignment
from sympy.codegen.ast import AddAugmentedAssignment
from sympy.printing.c import C99CodePrinter
# from sympy.codegen.cfunctions import *

output_dir = '/Users/zulianp/Desktop/code/utopia/utopia_fe/backend/kokkos/assembly/mech/generated'

from utopia_mech import *

mu, lmbda = symbols('mu lmbda')
dx = symbols('dx')


simplify_expressions = True

d = 3
# FE
grad_trial = trial_gradient(d)
grad_test = test_gradient(d)

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
# Strain energy function
###############################

Ogden = mu/2 *(I_C - d) - mu * log(J) + (lmbda/2) * (log(J))**2
Bower = mu/2 *(J**(-Rational(2, 3))*I_C - d) + (lmbda/2) * (J-1)**2
Wang  = mu/2 *(J**(-Rational(2, 3))*I_C - d) + (lmbda/2) * (J-1) # Does not work


mu_s = Rational(4,3)*mu
lmbda_s = lmbda + Rational(5,6)*lmbda
alpha = 1 + mu_s/lmbda_s - mu_s/(4*lmbda_s)

Smith = mu_s/2 * (I_C - d) + lmbda_s/2 * (J - alpha)**2 - mu_s/2 * log(I_C + 1)

W = simplify(Ogden)
# W = simplify(Bower)
# W = simplify(Wang)
# W = simplify(Smith)

P = first_piola(W, F)

# if False:
if True:

	#############################################
	# Energy
	#############################################
	energy = W * dx


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

		gradient_expression_list.append(AddAugmentedAssignment(symbols(f"lf[offset_i+{d1}]"), subsituted))

		for d2 in range(0, d):
			subsituted = tp.bilinear_subs("test", d1, "trial", d2, bilinear_form)

			if simplify_expressions:
				subsituted = simplify(subsituted)

			hessian_expression_list.append(AddAugmentedAssignment(symbols(f"bf[offset_ij+{d1*d + d2}]"), subsituted))
	

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


	generator.generate(hessian_expression_list, 'templates/utopia_tpl_elasticity_hessian.hpp', output_dir + '/generated_AutoHyperElasticityHessian.hpp')
	generator.generate(gradient_expression_list, 'templates/utopia_tpl_elasticity_gradient.hpp', output_dir + '/generated_AutoHyperElasticityGradient.hpp')
	generator.generate(energy_expression_list, 'templates/utopia_tpl_elasticity_energy.hpp', output_dir + '/generated_AutoHyperElasticityEnergy.hpp')

	generator.generate(full_expression_list, 'templates/utopia_tpl_elasticity.hpp', output_dir + '/generated_AutoHyperElasticity.hpp')

