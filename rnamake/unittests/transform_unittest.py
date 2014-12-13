import unittest
import rnamake.transform
import numpy as np
import numerical

class TransformUnittest(unittest.TestCase):

    def test_creation(self):
        t = rnamake.transform.Transform(np.eye(4))
        r = np.eye(3)
        d = np.array([0, 0, 0])
        if not numerical.are_matrices_equal(r, t.rotation()):
            self.fail("did not get rotation correct")

        if not numerical.are_points_equal(d, t.translation()):
            self.fail("did not get translation correct")

    def test_creation_2(self):
        r = np.eye(3)
        d = np.array([1, 1, 1])
        t = rnamake.transform.Transform(r, d)
        if not numerical.are_matrices_equal(r, t.rotation()):
            self.fail("did not save rotation correctly")
        if not numerical.are_points_equal(d, t.translation()):
            self.fail("did not save translation correctly")

def main():
    unittest.main()

if __name__ == '__main__':
    main()
