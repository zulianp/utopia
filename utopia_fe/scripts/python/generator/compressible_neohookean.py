# import sympy as sp
from sympy import *
from sympy.utilities.codegen import codegen
from sympy.utilities.codegen import CCodeGen
# from sympy.codegen.cfunctions import *

from utopia_mech import *


d = 2
# FE
grad_trial = trial_gradient(d)
# grad_trial = MatrixSymbol("grad_trial", 2, 2)

grad_test00, grad_test01, grad_test10, grad_test11 = symbols('grad_test_00 grad_test_01 grad_test_10 grad_test_11')
grad_test = Matrix(2,2,[grad_test00,grad_test01,grad_test10,grad_test11])
# grad_test = MatrixSymbol("grad_test", 2, 2)

# Nonlinear quantities
# 2D deformation gradient
f00, f01, f10, f11 = symbols('f00 f01 f10 f11')
C1, D1 = symbols('C1 D1')

F = deformation_gradient(d)
# F = MatrixSymbol("F", 2, 2)


J = det(F)

F_inv=Inverse(F)
F_inv_t = F_inv.T

I1 = trace(F)

print("Deformation Gradient")
print(F)

print("Det")
print(J)

print("I1")
print(I1)

W = simplify(C1 * (I1 - d - 2 * log(J)) + D1 * (J - 1)**2)

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
