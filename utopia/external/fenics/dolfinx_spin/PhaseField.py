from basix.ufl import element
from ufl import (Coefficient,FunctionSpace, Identity, TestFunction, TrialFunction,
                 FiniteElement, VectorElement, derivative, det, diff, dx, grad, ln,
                 tetrahedron, hexahedron, tr, variable, sym, inner, Mesh)

# Function spaces
disp_element = element("Lagrange", "tetrahedron", 1, shape=(3,))
phase_element = element("Lagrange", "tetrahedron", 1)

mesh = Mesh(disp_element)
C = FunctionSpace(mesh, phase_element)
V = FunctionSpace(mesh, disp_element)

# Trial and test functions
trial_u = TrialFunction(V)
test_u = TestFunction(V)

trial_c = TrialFunction(C)
test_c = TestFunction(C)

# Functions
u = Coefficient(V)
c = Coefficient(C)

# Kinematics
d = len(u)
I = Identity(d)

Gc = 1.7e-3
ls = 0.2
lambda_lame = 121.15
mu_lame = 80.77
kappa = lambda_lame+(2.*mu_lame/d)
deltaT = 5e-5
kdamage = 1e-9

def mc_bracket_p(v):
    return 0.5 * (v + abs(v))

def mc_bracket_n(v):
    return 0.5 * (v - abs(v))

def epsilon(u):
    return sym(grad(u))

def dev(epsilon):
    return (epsilon - ((1./d) * tr(epsilon)*I))

def g(c):
    return (1 - c)**2

def omega(c):
    return c**2

def c_omega(c):
    return 2

term1 = 0.5*kappa*(mc_bracket_p(tr(epsilon(u))))**2
dev_Epsilon = dev(epsilon(u))
term2 = mu_lame*inner(dev_Epsilon, dev_Epsilon)
elastic_energy_positive = ((1.-kdamage)* g(c) + kdamage)*(term1 + term2)
elastic_energy_negative = (0.5*kappa*(mc_bracket_n(tr(epsilon(u))))**2)

fracture_energy = (Gc/c_omega(c))*((omega(c)/ls) + ls*inner(grad(c), grad(c)))

# Displacement

disp_objective = ( elastic_energy_positive+elastic_energy_negative ) * dx
disp_gradient = derivative(disp_objective, u, test_u)
disp_hessian = derivative(disp_gradient, u, trial_u)

# Phase-field

phase_objective = fracture_energy * dx
phase_gradient = derivative(phase_objective, c, test_c)
phase_hessian = derivative(phase_gradient, c, trial_c)

# All outputs
forms = [disp_objective, disp_gradient, disp_hessian, phase_objective, phase_gradient, phase_hessian]
elements = [(disp_element), (phase_element)]
