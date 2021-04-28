import utopia as u 
import torch
import numpy as np


# Given a pytorch array, return an utopia vector
def pytorch_to_utopia(x):
	x_into_numpy = x.numpy()

	double_array = u.new_double_array(x_into_numpy.size)
	for i in range(0, x_into_numpy.size):
		u.double_array_setitem(double_array,i,x_into_numpy[i])

	comm = u.Communicator()
	size = x_into_numpy.size
	l = u.Layout(comm, size, size)

	u_vector = u.Vector()

	u_vector.serial_uconversion(double_array, size)

	u_vector.describe()
	return u_vector	
	# u.delete_double_array(a)


# Given a numpy array, return an utopia vector
def numpy_to_utopia(x):
	double_array = u.new_double_array(x.size)
	for i in range(0, x.size):
		u.double_array_setitem(double_array,i,x[i])

	comm = u.Communicator()
	size = x.size
	l = u.Layout(comm, size, size)

	u_vector = u.Vector()

	u_vector.serial_uconversion(double_array, size)

	u_vector.describe()
	return u_vector	
	# u.delete_double_array(a)