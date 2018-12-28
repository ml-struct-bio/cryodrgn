'''
Tools for dealing with SO(3) group and algebra
Adapted from https://github.com/pimdh/lie-vae
All functions are pytorch-ified
'''

import torch
from torch.distributions import Normal
import numpy as np

def map_to_lie_algebra(v):
    """Map a point in R^N to the tangent space at the identity, i.e.
    to the Lie Algebra
    Arg:
        v = vector in R^N, (..., 3) in our case
    Return:
        R = v converted to Lie Algebra element, (3,3) in our case"""

    # make sure this is a sample from R^3
    assert v.size()[-1] == 3

    R_x = v.new_tensor([[ 0., 0., 0.],
                        [ 0., 0.,-1.],
                        [ 0., 1., 0.]])

    R_y = v.new_tensor([[ 0., 0., 1.],
                        [ 0., 0., 0.],
                        [-1., 0., 0.]])

    R_z = v.new_tensor([[ 0.,-1., 0.],
                        [ 1., 0., 0.],
                        [ 0., 0., 0.]])

    R = R_x * v[..., 0, None, None] + \
        R_y * v[..., 1, None, None] + \
        R_z * v[..., 2, None, None]
    return R

def expmap(v):
    theta = v.norm(p=2, dim=-1, keepdim=True)
    # normalize K
    K = map_to_lie_algebra(v / theta)

    I = torch.eye(3, device=v.device, dtype=v.dtype)
    R = I + torch.sin(theta)[..., None]*K \
        + (1. - torch.cos(theta))[..., None]*(K@K)
    return R

def s2s1rodrigues(s2_el, s1_el):
    K = map_to_lie_algebra(s2_el)
    cos_theta = s1_el[...,0]
    sin_theta = s1_el[...,1]
    I = torch.eye(3, device=s2_el.device, dtype=s2_el.dtype)
    R = I + sin_theta[..., None, None]*K \
        + (1. - cos_theta)[..., None, None]*(K@K)
    return R

def s2s2_to_SO3(v1, v2):
    '''Normalize 2 3-vectors. Project second to orthogonal component.
    Take cross product for third. Stack to form SO matrix.'''
    u1 = v1
    e1 = u1 / u1.norm(p=2, dim=-1, keepdim=True).clamp(min=1E-5)
    u2 = v2 - (e1 * v2).sum(-1, keepdim=True) * e1
    e2 = u2 / u2.norm(p=2, dim=-1, keepdim=True).clamp(min=1E-5)
    e3 = torch.cross(e1, e2)
    return torch.stack([e1, e2, e3], 1)

def logsumexp(inputs, dim=None, keepdim=False):
    '''Numerically stable logsumexp.
    https://github.com/pytorch/pytorch/issues/2591

    Args:
        inputs: A Variable with any shape.
        dim: An integer.
        keepdim: A boolean.

    Returns:
        Equivalent of log(sum(exp(inputs), dim=dim, keepdim=keepdim)).
    '''
    # For a 1-D array x (any array along a single dimension),
    # log sum exp(x) = s + log sum exp(x - s)
    # with s = max(x) being a common choice.
    if dim is None:
        inputs = inputs.view(-1)
        dim = 0
    s, _ = torch.max(inputs, dim=dim, keepdim=True)
    outputs = s + (inputs - s).exp().sum(dim=dim, keepdim=True).log()
    if not keepdim:
        outputs = outputs.squeeze(dim)
    return outputs

def so3_entropy_old(w_eps, std, k=10):
    '''
    w_eps(Tensor of dim 3): sample from so3
    covar(Tensor of dim 3x3): covariance of distribution on so3
    k: 2k+1 samples for truncated summation
    '''
    # entropy of gaussian distribution on so3
    # see appendix C of https://arxiv.org/pdf/1807.04689.pdf
    theta = w_eps.norm(p=2)
    u = w_eps/theta # 3
    angles = 2*np.pi*torch.arange(-k,k+1,dtype=w_eps.dtype,device=w_eps.device) # 2k+1
    theta_hat = theta + angles # 2k+1
    x = u[None,:] * theta_hat[:,None] # 2k+1 , 3
    log_p = Normal(torch.zeros(3,device=w_eps.device),std).log_prob(x) # 2k+1,3
    clamp = 1e-3
    log_vol = torch.log((theta_hat**2).clamp(min=clamp)/(2-2*torch.cos(theta)).clamp(min=clamp)) # 2k+1
    log_p = log_p.sum(-1) + log_vol
    entropy = -logsumexp(log_p)
    return entropy

def so3_entropy(w_eps, std, k=10):
    '''
    w_eps(Tensor of dim Bx3): sample from so3
    std(Tensor of dim Bx3): std of distribution on so3
    k: Use 2k+1 samples for truncated summation
    '''
    # entropy of gaussian distribution on so3
    # see appendix C of https://arxiv.org/pdf/1807.04689.pdf
    theta = w_eps.norm(p=2, dim=-1, keepdim=True) # [B, 1]
    u = w_eps/theta # [B, 3]
    angles = 2*np.pi*torch.arange(-k,k+1,dtype=w_eps.dtype,device=w_eps.device) # 2k+1
    theta_hat = theta[:, None, :] + angles[:, None] # [B, 2k+1, 1]
    x = u[:,None,:] * theta_hat # [B, 2k+1 , 3]
    log_p = Normal(torch.zeros(3,device=w_eps.device),std).log_prob(x.permute([1,0,2])) # [2k+1, B, 3]
    log_p = log_p.permute([1,0,2]) # [B, 2k+1, 3]
    clamp = 1e-3
    log_vol = torch.log((theta_hat**2).clamp(min=clamp)/(2-2*torch.cos(theta_hat)).clamp(min=clamp)) # [B, 2k+1, 1]
    log_p = log_p.sum(-1) + log_vol.sum(-1) #[B, 2k+1]
    entropy = -logsumexp(log_p, -1)
    return entropy


