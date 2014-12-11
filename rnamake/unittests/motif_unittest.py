import unittest
import rnamake.motif


class MotifUnittest(unittest.TestCase):

    def test_creation(self):
        path = "/Users/josephyesselman/projects/REDESIGN/redesign/tests/p4p6"
        m = rnamake.motif.Motif(path)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
