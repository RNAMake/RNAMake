import numpy as np

class Transform(object):
    """
    Basic transform object, stores a 4x4 homogenous transform that contains
    both the rotation and translation of a transform

    .. code-block:: python

        #with full matrix
        >>> t = Transform(np.eye(4))
        >>> t.rotation()
        [[ 1.  0.  0.]
         [ 0.  1.  0.]
         [ 0.  0.  1.]]
        >>> print t.translation()
        [ 0.  0.  0.]

        # with rotation and translation
        >>> r = np.eye(3)
        >>> d = np.array([1, 1, 1])
        >>> t = Transform(r, d)
        >>> print t.rotation()
        [[ 1.  0.  0.]
         [ 0.  1.  0.]
         [ 0.  0.  1.]]
        >>> print t.translation()
        [ 1.  1.  1.]

    :attributes:

    `matrix` : np.array
        4x4 matrix that holds both the rotation and translation

    """
    __slots__ = ["matrix"]

    def __init__(self, *args):
        nargs = len(args)
        if nargs == 0:
            self.matrix = np.eye(4)
        elif nargs == 1:
            self.matrix = args[0]
        elif nargs == 2:
            self.matrix = np.eye(4)
            self.matrix[:3, :3] = args[0]
            self.matrix[:3, 3] = args[1]
        else:
            raise ValueError("invalid number of arguments")

    def rotation(self):
        """
        get rotation of this transform
        """
        return self.matrix[:3, :3]

    def translation(self):
        """
        get translation of this transform
        """
        return self.matrix[:3, 3]


