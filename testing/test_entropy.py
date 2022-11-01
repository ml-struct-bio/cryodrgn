import numpy as np
import torch

import cryodrgn.lie_tools

avg = []
std = torch.tensor([2.3407, 1.0999, 1.2962])
for _ in range(10):
    w_eps = torch.randn_like(std) * std
    e = cryodrgn.lie_tools.so3_entropy_old(w_eps, std)
    avg.append(e)
    print(e)

w_eps = torch.tensor([-0.46, -1.54, -1.96])
e = cryodrgn.lie_tools.so3_entropy_old(w_eps, std)
print(e)
avg.append(e)

# test new multi sample so3_entropy
w_eps = torch.tensor([-0.46, -1.54, -1.96])
w_eps2 = torch.tensor([1.0, 1.0, 1.0])
std = torch.tensor([2.3407, 1.0999, 1.2962])
std2 = torch.tensor([1.0, 1.0, 1.0])
e = cryodrgn.lie_tools.so3_entropy_old(w_eps2, std2)
print(e)
e = cryodrgn.lie_tools.so3_entropy(
    torch.stack([w_eps, w_eps2]), torch.stack([std, std2])
)
print(e)

a = torch.mean(torch.Tensor(avg))
print("average: {}".format(a))
print(np.log(8 * np.pi**2))
