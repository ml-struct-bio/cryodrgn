import torch


class Torch:

    array = torch.tensor

    def linspace(self, start, stop, num, **kwargs):
        endpoint = kwargs.pop('endpoint', True)
        if endpoint:
            return torch.linspace(start, stop, num, **kwargs)
        else:
            return torch.linspace(start, stop, num+1, **kwargs)[:-1]

    def __getattr__(self, item):
        """
        Catch-all method to to allow a straight pass-through \
        of any attribute that is not supported above.
        """

        return getattr(torch, item)
