import unittest
import random
import rnamake.vienna as vienna

def random_sequence():
    s = ""
    size = random.randrange(50, 200)
    e = "AUCG"
    for i in range(size):
        s += random.choice(e)
    return s


class ViennaUnittest(unittest.TestCase):

    def test_creation(self):
        v = vienna.Vienna()

    def test_fold(self):
        v = vienna.Vienna()
        results = v.fold(random_sequence())

    # TODO move to intergration
    def _test_fold_exhustive(self):
        v = vienna.Vienna()
        for i in range(100):
            results = v.fold(random_sequence())


def main():
    unittest.main()

if __name__ == '__main__':
    main()
