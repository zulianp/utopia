# import sympy as sp
from sympy import *
from sympy.utilities.codegen import codegen
from sympy.utilities.codegen import CCodeGen
# from sympy.codegen.cfunctions import *

from utopia_mech import *

mu, lmbda = symbols('mu lmbda')


d = 2
# FE
grad_trial = trial_gradient(d)
grad_test = test_gradient(d)

# Nonlinear quantities
F = deformation_gradient(d)
J = det(F)
F_inv = Inverse(F)
F_inv_t = F_inv.T
C = F.T*F
I_C = trace(C)

Ogden = mu/2 *(I_C - d) - mu * log(J) + (lmbda/2) * (log(J))**2
Bower = mu/2 *(I_C**(-Rational(2, 3)) - d) + (lmbda/2) * (J-1)**2
Wang  = mu/2 *(I_C**(-Rational(2, 3)) - d) + (lmbda/2) * (J-1)


mu_s = Rational(4,3)*mu
lmbda_s = lmbda + Rational(5,6)*lmbda
alpha = 1 + mu_s/lmbda_s - mu_s/(4*lmbda_s)

Smith = mu_s/2 * (I_C - d) + lmbda_s/2 * (J - alpha)**2 - mu_s/2 * log(I_C + 1)


# Strain energy function
# W = simplify(Ogden)
# W = simplify(Bower)
# W = simplify(Wang)
W = simplify(Smith)


P00 = diff(W, F[0, 0]);
P01 = diff(W, F[0, 1]);
P10 = diff(W, F[1, 0]);
P11 = diff(W, F[1, 1]);

P = Matrix(2,2,[P00,P01,P10,P11])

# if False:
if True:

	code_energy = ccode(W, standard='C89')
	print("Code energy:")
	print(code_energy)
	print("")

	print("Code linear form:")

	linear_form = simplify(trace(P.T * grad_test))
	code_lf =  ccode(linear_form, standard='C89')
	print(code_lf)


	print("Code gradient:")

	for i in range(0, d):
		for j in range(0, d):

			code_gradient =  ccode(P[i,j], standard='C89')
			print(f'{i},{j}) {code_gradient}')


	# for i in range(0, d):
	# 	for j in range(0, d):
	# 		for l in range(0, d):
	# 			for k in range(0, d):
	# 				Hij = diff(P[l, k], F[i, j]);

	# 				code_hessian =  ccode(Hij, standard='C89')
	# 				print(f'{i},{j},{l},{k}) {code_hessian}')

	print("")

	contraction = 0

	for i in range(0, d):
		for j in range(0, d):
	 		contraction += P[i,j] * grad_trial[i,j]


	print("Code stress")

	bilinear_form = 0

	for i in range(0, d):
		for j in range(0, d):
			Hij = simplify(diff(contraction, F[i, j]));

			bilinear_form += Hij * grad_test[i,j]

			Hijs = Hij;

			# Hijs = Hij.subs(grad_trial00, 0)
			# Hijs = Hijs.subs(grad_trial01, 1)
			# Hijs = Hijs.subs(grad_trial10, 0)
			# Hijs = Hijs.subs(grad_trial11, 0)
			# Hijs = simplify(Hijs)

			code_lin_stress = ccode(Hijs, standard='C89')
			print(f'{i},{j}) {code_lin_stress}')
			print("")

	print("")

	print("Code hessian")

	code_hessian = ccode(simplify(bilinear_form), standard='C89')
	print(f'{code_hessian}')
	print("")


# generator=CCodeGen(project='utopia', printer=None, preprocessor_statements=None, cse=True)
# r = generator.routine("neohookean", W, argument_sequence=None, global_vars=None)
# print(r)

# generator.dump_c((r,r), "neohookean.C", "neohookean", header=True, empty=True)

# [(c_name, c_code), (h_name, c_header)] = codegen(("neohookean", W), "C89", "test", header=False, empty=False)
# print(c_header)
# print(c_code)

# dr = diff(r*m_inv, m[0,1])
# # print(dr)


# dr=simplify(dr)


# print("C---------\n")
# print(cc)

# print("C---------\n")

# # https: // docs.sympy.org/latest/modules/diffgeom.html
