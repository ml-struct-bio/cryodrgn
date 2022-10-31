import torch

from cryodrgn import lie_tools


def test_so3_entropy():
    entropy = lie_tools.so3_entropy(
        w_eps=torch.stack(
            [torch.tensor([-0.46, -1.54, -1.96]), torch.tensor([1.0, 1.0, 1.0])]
        ),
        std=torch.stack(
            [torch.tensor([2.3407, 1.0999, 1.2962]), torch.tensor([1.0, 1.0, 1.0])]
        ),
    )
    assert torch.allclose(entropy, torch.tensor([5.3786, 3.9993]))
