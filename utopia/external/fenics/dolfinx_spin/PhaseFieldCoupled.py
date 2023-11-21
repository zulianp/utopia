from ufl import (Coefficient, Identity, 
                 TestFunction, TrialFunction, TestFunctions, TrialFunctions,
                 FiniteElement, VectorElement, derivative, det, diff, dx, grad, ln,
                 tetrahedron, tr, variable, sym, inner, MixedElement, split)

vector_element = VectorElement("Lagrange", tetrahedron, 1)
element = FiniteElement("Lagrange", tetrahedron, 1)
mixed = MixedElement([vector_element, element])

trial = TrialFunction(mixed)
test = TestFunction(mixed)

x = Coefficient(mixed)
u, c = split(x)

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

objective = (fracture_energy + elastic_energy_positive + elastic_energy_negative ) * dx
gradient = derivative(objective, x, test)
hessian = derivative(gradient, x, trial)

# All outputs

forms = [objective, gradient, hessian]
elements = [(mixed)]
