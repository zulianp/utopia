from basix.ufl import element
from ufl import (Coefficient, FunctionSpace, Identity, Mesh, TestFunction,
                 TrialFunction, derivative, det, diff, dx, grad, ln, tr,
                 variable)


e = element("Lagrange", "tetrahedron", 1, shape=(3,))
mesh = Mesh(e)
V = FunctionSpace(mesh, e)


du = TrialFunction(V)
v = TestFunction(V)

u = Coefficient(V)
u_old = Coefficient(V)

d = len(u)
I = Identity(d)
F = variable(I + grad(u))
C = F.T * F

F_old = variable(I + grad(u_old))
C_old = F_old.T * F_old

Ic = tr(C)
J = det(F)

dt = 0.01
E = 10.0
nu = 0.3
mu = E / (2 * (1 + nu))
lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
eta = 0.05

elastic_energy = (mu / 2) * (Ic - 3) - mu * ln(J) + (lmbda / 2) * (ln(J))**2
viscous_energy = eta/dt  * (tr(C*C) - 2 * tr(F * C_old * F.T))

energy = (elastic_energy + viscous_energy) * dx
gradient = derivative(energy, u, v)
hessian = derivative(gradient, u, du)

forms = [energy, gradient, hessian]
elements = [e]
expressions = []
