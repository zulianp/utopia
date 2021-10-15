import utopia as u
import torch

x = torch.rand(5, 3)
print(x)


u.init()

b = u.Vector()

b.set(5,4)

A = u.SparseMatrix()

A.print_info()

u.finalise()

#cd utopia/utopia/build/scripting && python3 testing_swig.py