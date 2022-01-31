
# https://www.youtube.com/watch?v=5jzIVp6bTy0&ab_channel=Enthought

# https://docs.sympy.org/latest/modules/solvers/pde.html

# https://cardona.co/math/2018/09/16/geodesics.html

# import sympy as sp
from sympy import *
from sympy.utilities.codegen import codegen
from sympy.utilities.codegen import CCodeGen
from sympy.printing.c import C99CodePrinter
from sympy.codegen.ast import Assignment


# class UtopiaPrinter(C99CodePrinter):
# 	def _print_list(self, list_of_exprs):
# 		if(all(isinstance(x, ImmutableMatrix) for x in list_of_exprs)):
# 			sub_exprs, simplified_exprs = cse(list_of_exprs)
# 			lines = []
# 			for var, sub_expr in sub_exprs:
# 				# ass = Assignment(var, sub_expr.xreplace(state_array_map))
# 				ass = Assignment(var, sub_expr) # replace?
# 				lines.append('ValueType ' + self._print(ass))

# 			for mat in simplified_exprs:
# 			# 	# lines.append(self._print(mat.xreplace(state_array_map)))
# 				print(mat)
# 				# g_result =MatrixSymbol("uh", *mat.shape)
# 				# lines.append(self._print(mat,assign_to=g_result))
# 			return '\n'.join(lines)
# 		else:
# 			return super()._print_list(list_of_exprs)

# 	def _print_ImmutableDenseMatrix(self, expr):
# 		if expr.shape[1] > 1:
# 			M = MatrixSymbol("matrix_result", *expr.shape)
# 		else:
# 			M = Matrix("vector_result", *expr.shape)
# 		return self._print(Assignment(M, expr))


	# def _print_ImmutableDenseMatrix(self, expr):
	# 	sub_exprs, simplified = cse(expr)
	# 	lines = []
	# 	for var, sub_expr in sub_exprs:
	# 		lines.append('ValueType ' + self._print(Assignment(var, sub_expr)))
	# 	M = MatrixSymbol('M', *expr.shape)
	# 	return '\n'.join(lines) + '\n' + self._print(Assignment(M, expr))


class ReferenceTriangle:
	def f0(self, x):
		return 1-x[0]-x[1]
	def f1(self, x):
		return x[0]
	def f2(self, x):
		return x[1]
	def f(self, x):
		return Matrix([self.f0(x), self.f1(x), self.f2(x)])


class Triangle:
	def __init__(self,p0,p1,p2):
		self.b = p0
		self.A = Matrix([[p1[0]-p0[0],
						 p2[0]-p0[0]],
						[p1[1]-p0[1],
						 p2[1]-p0[1]]])

		self.f0s = Function('f0');
		self.f1s = Function('f1');
		self.f2s = Function('f2');

	def f0(self, x):
		return self.f0s(x[0], x[1])
	def f1(self, x):
		return self.f1s(x[0], x[1])
	def f2(self, x):
		return self.f2s(x[0], x[1])
	def f(self, x):
		return Matrix([self.f0(x), self.f1(x), self.f2(x)])


	# def f0(self, x):
	# 	x_hat = self.inverse_transform(x)
	# 	return ReferenceTriangle().f0(x_hat)

	# def f1(self, x):
	# 	x_hat = self.inverse_transform(x)
	# 	return ReferenceTriangle().f1(x_hat)

	# def f2(self, x):
	# 	x_hat = self.inverse_transform(x)
	# 	return ReferenceTriangle().f2(x_hat)

	# def f(self, x):
	# 	x_hat = self.inverse_transform(x)
	# 	return ReferenceTriangle().f(x_hat)

	# def inverse_transform(self, x):
	# 	x_hat = Inverse(self.A) * (x - self.b)
	# 	return x_hat

	# def transform(self, x_hat):
	# 	x = self.A * x_hat + self.b
	# 	return x




init_printing()

d = 2


x,y,t,c = symbols('x y t c')
# alpha = symbols('alpha')
# D = MatrixSymbol("D", 2, 2)
# A = Matrix([[1, 2], [3, 4]])

# J = Matrix([sin(x)*cos(t), sin(x)]).jacobian([x, t])

# print(J)
# 
# f = Function('f') # Parameters can be passed later
# u = Function('u')(x, t)
# g = f(x) + x *alpha * Rational(1, 2)

# print(g)

# print(diff(g, x))

# eq=Eq(diff(u, t, t), c**2*diff(u, x, 2))

# print(eq)
p0x, p1x, p2x = symbols('p0x p1x p2x')
p0y, p1y, p2y = symbols('p0y p1y p2y')

p0 = Matrix([p0x, p0y])
p1 = Matrix([p1x, p1y])
p2 = Matrix([p2x, p2y])


# fe = Triangle(Matrix([0,0]),Matrix([2,0]),Matrix([0, 3]))
fe = Triangle(p0, p1, p2);
# v = simplify(diff(fe.f0(Matrix([x,t])), Matrix([x,t])))
# v = simplify(diff(fe.f(Matrix([x,y])), x))

# x_hat,y_hat = symbols('x_hat y_hat')
# p_hat = Matrix([x_hat,y_hat])
# p = fe.transform(p_hat)

# print(p)

p = Matrix([x,y])
v = fe.f(p).jacobian(p)


# print(v)
# print(f"Code: {v.shape}")

# Each basis function
# for i in range(0, 3):
# 		print(f'df(p{i})_dx = {ccode(v[i,0])}')
# 		print(f'df(p{i})_dy = {ccode(v[i,1])}')

u0, u1, u2 =  symbols('u0 u1 u2')

uh = simplify(u0 * fe.f0(p) + u1 * fe.f1(p) + u2 * fe.f2(p))

# print(f'uh = {uh}')

grad_uh = diff(simplify(uh),p)

grad_uh = grad_uh.subs(Derivative(fe.f0(p), x), symbols('g0x'))
grad_uh = grad_uh.subs(Derivative(fe.f0(p), y), symbols('g0y'))

grad_uh = grad_uh.subs(Derivative(fe.f1(p), x), symbols('g1x'))
grad_uh = grad_uh.subs(Derivative(fe.f1(p), y), symbols('g1y'))

grad_uh = grad_uh.subs(Derivative(fe.f2(p), x), symbols('g2x'))
grad_uh = grad_uh.subs(Derivative(fe.f2(p), y), symbols('g2y'))

uh = uh.subs(fe.f0(p), symbols('f0'))
uh = uh.subs(fe.f1(p), symbols('f1'))
uh = uh.subs(fe.f2(p), symbols('f2'))

print(grad_uh)



# print(f'duh_dx = {ccode(grad_uh[0])}')
# print(f'duh_dy = {ccode(grad_uh[1])}')

u_result = symbols('uh')
g_result =MatrixSymbol("guh", *grad_uh.shape)
sub_expr, simpl_expr = cse((grad_uh, uh))


# help(cse)

# print("temp")

lines = []
printer = C99CodePrinter()

for var,expr in sub_expr:
	lines.append(f'T {var} = {printer.doprint(expr)};')

# print(simpl_expr[0][0])
# print(simpl_expr[0][1])

# print("store")



for v in simpl_expr:
	if isinstance(v, ImmutableMatrix):
		lines.append(printer.doprint(v, assign_to=g_result))
	else:
		lines.append(printer.doprint(v, assign_to=u_result))



code_string='\n'.join(lines)
code_tpl = """\
template<typename T>
void f(const T *u, const T *phi, const T*grad_phi[3], constT &uh, T * guh) 
{{ 

T u0 = u[0];
T u1 = u[1];
T u2 = u[2];

T f0 = phi[0];
T f1 = phi[1];
T f2 = phi[2];

T g0x = grad_phi[0][0];
T g1x = grad_phi[1][0];
T g2x = grad_phi[2][0];

T g0y = grad_phi[0][1];
T g1y = grad_phi[1][1];
T g2y = grad_phi[2][1];

{code}
}}
"""

fun = code_tpl.format(code=code_string)
print(fun)


# with open('f.c', 'w') as f:
# 	f.write(fun)
 

# printer = C99CodePrinter()
# code=printer.doprint(grad_uh, assign_to=g_result)
# printer = UtopiaPrinter()
# code = printer.doprint([uh,grad_uh])
# 
# print(code)
# l = latex(grad_uh)
# print(l)

