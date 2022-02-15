# import sympy as sp
from sympy import *
import os
import sys
import rich

from utopia_mech import *

console = rich.get_console()

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
class SaintVenantKirchoff(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		mu, lmbda = symbols('mu lambda')
		self.params = [(mu, 1.0), (lmbda, 1.0)]

		C = self.C
		J = self.J
		E = (C - eye(d, d))/2

		self.fun = lmbda/2 * trace(E)**2 + mu * trace(E*E)
		self.name = 'SaintVenantKirchoff'


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

# TEST (https://shaddenlab.gitlab.io/fenicsmechanics/chapters/demos-all.html)
# https://www2.karlin.mff.cuni.cz/~hron/fenics-tutorial/elasticity/doc.html
# https://abaqus-docs.mit.edu/2017/English/SIMACAEMATRefMap/simamat-c-hyperelastic.htm

class IncompressibleMooneyRivlin(IncompressibleHyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		I_C = self.I_C
		J = self.J
		C = self.C
		p = self.p

		C1, C2, a = symbols('C1 C2 a')
		self.params = [(C1, 1.0), (C2, 1.0), (a, 1.)]

		CikxCki = 0
		for i in range(0, d):
			for k in range(0, d):
				CikxCki += C[i,k] * C[k,i]

		I_C2 = Rational(1, 2) * (I_C**2 - CikxCki)

		# self.fun = C1 * (I_C - d) + C2 * (I_C2 - d) + p/2 * (J - 1)**2
		# self.fun = C1 * (I_C - d) + C2 * (I_C2 - d) + p * (J - 1) - 0.5 * a * p * p
		self.fun = C1 * (I_C - d) + C2 * (I_C2 - d) + p * (J - 1)
		self.name = 'IncompressibleMooneyRivlin'

def generate_materials(d,simplify_expressions):
	output_dir = f'../../../backend/kokkos/assembly/mech/generated/{d}D'
	# models = [NeoHookeanOgden(d), NeoHookeanBower(d), NeoHookeanWang(d), NeoHookeanSmith(d), Fung(d), MooneyRivlin(d), SaintVenantKirchoff(d), IncompressibleMooneyRivlin(d)]
	# models = [ NeoHookeanSmith(d), Fung(d)]
	# models = [NeoHookeanOgden(d)]
	# models = [Fung(d)]
	models = [IncompressibleMooneyRivlin(d)] 
	# models = [SaintVenantKirchoff(d)]

	for m in models:
		m.generate_files(output_dir, simplify_expressions)

def main(args):
	generate_materials(2,True)
	# generate_materials(3, False)
	# generate_materials(3, True)

if __name__ == '__main__':
    main(sys.argv[1:])
