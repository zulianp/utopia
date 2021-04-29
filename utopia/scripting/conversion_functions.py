import utopia as u 
import torch
import numpy as np


# Given a pytorch mono-dimensional tensor, return an utopia vector
def pytorch_to_utopia(x):
	x_into_numpy = x.detach().numpy().flatten()

	double_array = u.new_double_array(x_into_numpy.size)
	for i in range(0, x_into_numpy.size):
		u.double_array_setitem(double_array,i,np.double(x_into_numpy[i]))

	u_vector = u.Vector()
	size = x_into_numpy.size

	u_vector.serial_uconversion(double_array, size)
	u.delete_double_array(double_array)

	# u_vector.describe()
	return u_vector	
	


# Given a numpy array, return an utopia vector
def numpy_to_utopia(x):
	double_array = u.new_double_array(x.size)
	for i in range(0, x.size):
		u.double_array_setitem(double_array,i,x[i])

	u_vector = u.Vector()
	size = x.size

	u_vector.serial_uconversion(double_array, size)
	u.delete_double_array(double_array)

	u_vector.describe()
	return u_vector	


# Given a utopia vector, return a pytorch mono-dimensional tensor 
def utopia_to_pytorch(vec, size):
	values = np.zeros((size,))
	
	for i in range(size):
		values[i] = vec.get(i)
	output = torch.from_numpy(values)

	return output













