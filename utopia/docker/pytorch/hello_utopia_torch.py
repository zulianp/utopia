import torch
import utopia
import numpy

utopia.init()

print("torch version: ")
print(torch.version.git_version)


a = numpy.array([1, 2, 3])
utopia.print_array(a)

utopia.finalize()