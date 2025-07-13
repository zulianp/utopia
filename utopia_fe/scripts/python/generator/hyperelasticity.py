# import sympy as sp
from sympy import *
import os
import sys
import rich

from utopia_mech import *


console = rich.get_console()


# class Kulkarni(HyperElasticModel):
# 	def __init__(self, d):
# 		super().__init__(d)

# 		C = self.C
# 		I1 = self.I_C
# 		I2 = se.rational(1, 2) * (trace(C)**2 - trace(C*C))
# 		I3 = det(C)
# 		I4 = 1

# 		F = self.F

# 		self.L = velocity_gradient(d)
# 		D = se.rational(1, 2) * (L + L.T)

# 		C_dot = 2 * F.T * D * F
# 		C_dot_symb = symbolic_matrix(d, 'C_dot')

# 		J2 = se.rational(1, 2) * trace(C_dot_symb**2)
# 		J5 = J2 # Isotropic

# 		Wv = se.rational(1, 2) * mu2 * J2 * (I1 - 3)**n1 + mu3*J5 * (K4 - 7) ** n2


		# J = self.J



		# mu, mu1, lmbda = se.symbols('mu mu1 lambda')
		# self.params = [(mu, 1.0), (lmbda,1.0)]


		# K4 = I2 + 2 * I1*I4 - 2*I5
		# self.fun =  se.rational(mu, 2) * (I1 - 3) + mu1 * (K4 - 7)**q


		# self.fun = mu/2 *(I_C - d) - mu * log(J) + (lmbda/2) * (log(J))**2
		# self.name = 'Kulkarni'
		# self.use_default_parameter_reader = True

# 1D equation (39)
# 3D eqaution (50, 51)
# History variable h is just evaluated for the update (explicit integrtaion)
# 
# class MaxwellElement(HyperElasticModel):
# 	def __init__(self, d):
# 		super().__init__(d)

# 		I_C = self.I_C
# 		J = self.J

# 		mu, lmbda = se.symbols('mu lambda')
# 		self.params = [(mu, 1.0), (lmbda,1.0)]

# 		self.fun = mu/2 *(I_C - d) - mu * log(J) + (lmbda/2) * (log(J))**2
# 		self.name = 'MaxwellElement'
# 		self.use_default_parameter_reader = True

###############################
# NeoHookean models
###############################

class NeoHookeanOgden(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		I_C = self.I_C
		J = self.J

		mu, lmbda = se.symbols('mu lambda')
		self.params = [(mu, 1.0), (lmbda,1.0)]

		self.fun = mu/2 *(I_C - d) - mu * log(J) + (lmbda/2) * (log(J))**2
		self.name = 'NeoHookeanOgden'
		self.use_default_parameter_reader = True

# https://onlinelibrary.wiley.com/doi/full/10.1002/cnm.2945
class NeoHookeanSiguenza(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		I_C = self.I_C
		J = self.J

		I_mod = J**(-2/3) * I_C

		G, K = se.symbols('G K')
		self.params = [(G, 1.0), (K, 1.0)]

		self.fun = G/2 *(I_mod - d) + (K/2) * (log(J))**2
		self.name = 'NeoHookeanSiguenza'
		self.use_default_parameter_reader = False

class NeoHookeanBower(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		I_C = self.I_C
		J = self.J

		mu, lmbda = se.symbols('mu lambda')
		self.params = [(mu, 1.0), (lmbda, 1.0)]

		self.fun = mu/2 *(J**(-se.rational(2, 3))*I_C - d) + (lmbda/2) * (J-1)**2
		self.name = 'NeoHookeanBower'
		self.use_default_parameter_reader = True

class NeoHookeanWang(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		I_C = self.I_C
		J = self.J

		mu, lmbda = se.symbols('mu lambda')
		self.params = [(mu, 1.0), (lmbda, 1.0)]

		self.fun = mu/2 *(J**(-se.rational(2, 3))*I_C - d) + (lmbda/2) * (J-1) 
		self.name = 'NeoHookeanWang'
		self.use_default_parameter_reader = True

class NeoHookeanSmith(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		I_C = self.I_C
		J = self.J

		mu, lmbda = se.symbols('mu lambda')
		self.params = [(mu, 1.0), (lmbda, 1.0)]

		mu_s = se.rational(4,3)*mu
		lmbda_s = lmbda + se.rational(5,6)*lmbda
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

		a, b, c, k = se.symbols('a b c k')
		self.params = [(a, 1.0), (b, 1.0), (c, 1.0), (k, 1)]

		C = self.C
		J = self.J
		E = (C - eye(d, d))/2

		# Isotropic
		trE = trace(E)
		q = a * trE
		Q = b * trE

		self.fun = se.rational(1,2) * (q + c * (exp(Q) - 1)) + (k/2)*(J - 1)**2
		self.name = 'Fung'

###############################
class SaintVenantKirchoff(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		mu, lmbda = se.symbols('mu lambda')
		self.params = [(mu, 1.0), (lmbda, 1.0)]

		C = self.C
		J = self.J
		E = (C - eye(d, d))/2

		self.fun = lmbda/2 * trace(E)**2 + mu * trace(E*E)
		self.name = 'SaintVenantKirchoff'

class Yeoh(HyperElasticModel):
	def __init__(self, d, n):
		super().__init__(d)

		C0 = []
		C1 = []
		self.params = []

		for i in range(0, n):
			C0.append(se.symbols(f"C0_{i}"))
			C1.append(se.symbols(f"C1_{i}"))

			self.params.append((C0[i], 1.))
			self.params.append((C1[i], 1.))

		C = self.C
		J = self.J

		I1 = J**se.rational(-2,3) * trace(C)

		self.fun = 0

		for i in range(0, n):
			self.fun += C0[i] * (I1 - d)**(i+1) +  C1[i] * (J - 1)**(2 * (i+1))

		if n == 2:
			self.name = "Yeoh"
		else:
			self.name = f"Yeoh{n}"

###############################
# Mooney-Rivlin (https://en.wikipedia.org/wiki/Mooney%E2%80%93Rivlin_solid)
###############################
class MooneyRivlin(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		C1, C2, K = se.symbols('C1 C2 K')
		self.params = [(C1, 0.083), (C2, 0.083), (K, 166.67)]

		J = self.J
		C = self.F * self.F.T # Left Cauchy-Green

		I1 = se.trace(C)
		I2 = se.rational(1, 2) * (I1 ** 2 - se.trace(C*C))

		Jm23 = J**se.rational(-2, 3)
		Jm43 = J**se.rational(-4, 3)
		I1_bar = Jm23 * I1
		I2_bar = Jm43 * I2

		I1_ref = d
		I2_ref = d

		if d == 2:
			I2_ref = 1

		W_iso =  C1 * (I1_bar - I1_ref) + C2 * (I2_bar - I2_ref)
		W_vol = (K/2) * log(J)**2
		self.fun = W_iso + W_vol
		self.name = 'MooneyRivlin'

class IncompressibleMooneyRivlin(HyperElasticModel):
	def __init__(self, d):
		super().__init__(d)

		C1, C2, K = se.symbols('C1 C2 K')
		self.params = [(C1, 0.083), (C2, 0.083), (K, 166.67)]

		J = self.J
		C = self.F * self.F.T # Left Cauchy-Green

		I1 = se.trace(C)
		I2 = se.rational(1, 2) * (I1 ** 2 - se.trace(C*C))

		Jm23 = J**se.rational(-2, 3)
		Jm43 = J**se.rational(-4, 3)
		I1_bar = Jm23 * I1
		I2_bar = Jm43 * I2

		I1_ref = d
		I2_ref = d

		if d == 2:
			I2_ref = 1

		W_iso =  C1 * (I1_bar - I1_ref) + C2 * (I2_bar - I2_ref)
		W_vol = (K/2) * (J - 1)**2
		self.fun = W_iso + W_vol
		self.name = 'IncompressibleMooneyRivlin'

class GuccioneCosta(HyperElasticModel):
	def __init__(self):
		super().__init__(3)

		I_C = self.I_C
		J = self.J

		Id = se.eye(3)
		E = se.rational(1, 2) * (self.C - Id)

		# mu sometimes is called C
		mu, bf, bt, bfs, k = se.symbols('mu b_f b_t b_fs k')
		self.params = [(mu, 2000), (bf, 8), (bt, 2), (bfs, 4), (k, 10)]

		Q = bf   *  E[0, 0]**2
		Q += bt  * (E[1, 1]**2 + E[2, 2]**2 + E[1, 2]**2 + E[2, 1]**2) 
		Q += bfs * (E[0, 1]**2 + E[1, 0]**2 + E[0, 2]**2 + E[2, 0]**2)

		self.fun = mu/2 * (se.exp(Q) - 1) + k * (J - 1)**2
		self.name = 'GuccioneCosta'
		self.use_default_parameter_reader = False

# TEST (https://shaddenlab.gitlab.io/fenicsmechanics/chapters/demos-all.html)
# https://www2.karlin.mff.cuni.cz/~hron/fenics-tutorial/elasticity/doc.html
# https://abaqus-docs.mit.edu/2017/English/SIMACAEMATRefMap/simamat-c-hyperelastic.htm

# class IncompressibleMooneyRivlin(IncompressibleHyperElasticModel):
# 	def __init__(self, d):
# 		super().__init__(d)

# 		I_C = self.I_C
# 		J = self.J
# 		C = self.C
# 		p = self.p

# 		C1, C2, a = se.symbols('C1 C2 a')
# 		self.params = [(C1, 1.0), (C2, 1.0), (a, 1.)]

# 		CikxCki = 0
# 		for i in range(0, d):
# 			for k in range(0, d):
# 				CikxCki += C[i,k] * C[k,i]

# 		I_C2 = se.rational(1, 2) * (I_C**2 - CikxCki)

# 		# self.fun = C1 * (I_C - d) + C2 * (I_C2 - d) + p/2 * (J - 1)**2
# 		# self.fun = C1 * (I_C - d) + C2 * (I_C2 - d) + p * (J - 1) - 0.5 * a * p * p
# 		self.fun = C1 * (I_C - d) + C2 * (I_C2 - d) + p * (J - 1)
# 		self.name = 'IncompressibleMooneyRivlin'

def generate_materials(d,simplify_expressions):
	output_dir = f'../../../backend/kokkos/assembly/mech/generated/{d}D'
	# output_dir = f'./workspace/{d}D'
	models = [NeoHookeanOgden(d), NeoHookeanBower(d), NeoHookeanWang(d), NeoHookeanSmith(d), Fung(d), MooneyRivlin(d), SaintVenantKirchoff(d), Yeoh(d, 2),]
	# models = [MooneyRivlin(d)] 
	models = [IncompressibleMooneyRivlin(d)] 

	# models.append(Yeoh(d, 3))
	# models = [ NeoHookeanSmith(d), Fung(d)]
	# models = [NeoHookeanOgden(d)]
	# models = [Fung(d)]
	# models = [IncompressibleMooneyRivlin(d)] 
	# models = [Yeoh(d, 3)]
	 
	# models = [Yeoh(d, 1)]
	# models = [SaintVenantKirchoff(d)]

	for m in models:
		m.generate_files(output_dir, simplify_expressions)

def main(args):
	# generate_materials(2,True)
	generate_materials(2, False)
	generate_materials(3, False)
	# generate_materials(3, True)

	output_dir = f'../../../backend/kokkos/assembly/mech/generated/3D'
	# GuccioneCosta().generate_files(output_dir, False)
	# NeoHookeanSiguenza(2).generate_files(output_dir, False)
	# NeoHookeanSiguenza(3).generate_files(output_dir, False)
	# MooneyRivlin(3).generate_files(output_dir, False)
	# MooneyRivlin(2).generate_files(output_dir, False)

if __name__ == '__main__':
    main(sys.argv[1:])
