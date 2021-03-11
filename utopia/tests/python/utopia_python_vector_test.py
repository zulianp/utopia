import utopia as u 
import torch
import unittest

class TestUtopiaVector(unittest.TestCase):

	def test_vector_cration(self):
		u.init()
		comm = u.Communicator()
		l = u.Layout(comm, 10, 10)
		a = u.Vector()
		a.values(l, 10)
		b = u.Vector()
		b.values(l, 10)
		a.equals(b,0)



    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_isupper(self):
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)

if __name__ == '__main__':
    unittest.main()
