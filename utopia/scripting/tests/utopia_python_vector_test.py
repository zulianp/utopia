
import utopia as u 
import torch
import unittest

class TestUtopiaVector(unittest.TestCase):

	def test_vector_creation(self):
		u.init()
		
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
		u.init()
		
		comm = u.Communicator()
		l = u.Layout(comm, 10, 10)
		
		a = u.Vector()
		a.values(l, 10)
		
		b = u.Vector()
		b.values(l, 10)

		a.axpy(0.1, b)

		test = u.Vector()
		test.values(l,11)

		self.assertTrue(a.equals(test,0), True)

		c = u.Vector()
		c.values(l,10)

		d = u.Vector()
		d.values(l,5)

		c.axpy(0.1, d)

		self.assertFalse(c.equals(test,0), False)

	def test_vector_dot_product(self):
		u.init()
		
		comm = u.Communicator()
		l = u.Layout(comm, 10, 10)
		
		a = u.Vector()
		a.values(l, 10)	

		b = u.Vector()
		b.values(l, 10)	

		self.assertTrue(a.dot(b), 1000.0)

if __name__ == '__main__':
    unittest.main()