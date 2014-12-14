import unittest
import rnamake.base

class BaseUnittest(unittest.TestCase):

    def test_creation(self):
        try:
            base = rnamake.base.Base()
        except:
            self.fail("got an error I did not expect")

    def _test_option(self):
        base = rnamake.base.Base()
        base.option("test", 1)

        if base.option("test") != 1:
            self.fail("did not get option back")


def main():
    unittest.main()

if __name__ == '__main__':
    main()
