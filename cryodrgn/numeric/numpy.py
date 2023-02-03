import numpy as np


class Numpy:

    asnumpy = staticmethod(lambda x: x)

    def __getattr__(self, item):
        """
        Catch-all method to to allow a straight pass-through \
        of any attribute that is not supported above.
        """

        return getattr(np, item)
