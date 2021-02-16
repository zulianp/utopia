import utopia as u
import torch

x = torch.rand(5, 3)
print(x)


u.init()

b = u.Vector()
b.set(1,10)
b.set(3,10)
b.disp()

A = u.SparseMatrix()

A.print_info()

u.finalise()

#cd utopia/utopia/build/scripting && python3 testing_swig.py
#import utopia.build.scripting.utopia