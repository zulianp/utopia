
import utopia as u 
from conversion_functions import *
import torch
import unittest
import numpy as np

class TestUtopiaVector(unittest.TestCase):

	def test_vector_creation(self):
		
		comm = u.Communicator()
		l = u.Layout(comm, 10, 10)
		
		a = u.Vector()
		a.values(l, 10)
		
		b = u.Vector()
		b.values(l, 10)
		
		self.assertTrue(a.equals(b,0), True)

		c = u.Vector()
		c.values(l,10)

		d = u.Vector()
		d.values(l,5)

		self.assertFalse(c.equals(d,0), False)

	def test_vector_addition(self):
		
		comm = u.Communicator()
		l = u.Layout(comm, 10, 10)
		
		a = u.Vector()
		a.values(l, 10)
		
		b = u.Vector()
		b.values(l, 10)

		a.axpy(0.1, b)

		test = u.Vector()
		test.values(l,11)

		c = u.Vector()
		c.values(l,10)

		d = u.Vector()
		d.values(l,5)

		c.axpy(0.1, d)

		self.assertTrue(a.equals(test,0), True)
		self.assertFalse(c.equals(test,0), False)

	def test_vector_dot_product(self):
		
		comm = u.Communicator()
		l = u.Layout(comm, 10, 10)
		
		a = u.Vector()
		a.values(l, 10)	

		b = u.Vector()
		b.values(l, 10)	

		self.assertTrue(a.dot(b), 1000.0)

	def test_vector_scale(self):
		comm = u.Communicator()
		l = u.Layout(comm, 10, 10)
		
		a = u.Vector()
		a.values(l, 10)

		b = u.Vector()
		b.values(l, 20)	

		c = u.Vector()
		c.values(l, 3)	

		a.scale(2)

		self.assertTrue(a.equals(b,0), True)
		self.assertFalse(a.equals(c,0), False)


	def test_numpy_to_utopia(self):
		comm = u.Communicator()
		l = u.Layout(comm, 10, 10)
		
		a = u.Vector()
		a.values(l, 10)

		c = u.Vector()
		c.values(l, 11)

		b = np.full((10, 1), 10.0)

		b_ut = numpy_to_utopia(b)

		self.assertTrue(a.equals(b_ut,0), True)
		self.assertFalse(c.equals(b_ut,0), False)


	def test_pytorch_to_utopia(self):
		comm = u.Communicator()
		l = u.Layout(comm, 10, 10)
		
		a = u.Vector()
		a.values(l, 10)

		c = u.Vector()
		c.values(l, 11)

		b = np.full((10, 1), 10.0)
		b_torch = torch.from_numpy(b)

		b_ut = pytorch_to_utopia(b_torch)

		self.assertTrue(a.equals(b_ut,0), True)
		self.assertFalse(c.equals(b_ut,0), False)	


	def test_pytorch_to_utopia(self):
		comm = u.Communicator()
		l = u.Layout(comm, 10, 10)
		
		a = u.Vector()
		a.values(l, 10)

		a_np =  np.full((10, 1), 10)
		a_torch = torch.from_numpy(a_np)

		utopia_to_pytorch(a, a_torch)

		b = np.full((10, 1), 10)
		b_torch = torch.from_numpy(b)

		c = np.full((10, 1), 10)
		c_torch = torch.from_numpy(c)

		self.assertTrue(torch.equal(a_torch, b_torch), True)
		self.assertTrue(torch.equal(a_torch, c_torch), False)
	


if __name__ == '__main__':
	u.init()
	unittest.main()
	u.finalize()

























