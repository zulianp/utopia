import torch
import utopia

utopia.init()

print("torch version: ")
print(torch.version.git_version)

utopia.finalize()