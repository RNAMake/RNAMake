import unittest
import rnamake.basic_io
import numerical

class BasicIoUnittest(unittest.TestCase):

    def test_point_to_str(self):
        p = [1, 0, 1]
        s = rnamake.basic_io.point_to_str(p)
        p1 = rnamake.basic_io.str_to_point(s)

        if not numerical.are_points_equal(p, p1):
            self.fail("points are not the same")


def main():
    unittest.main()

if __name__ == '__main__':
    main()
