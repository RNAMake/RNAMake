import unittest
import rnamake.base

class BaseUnittest(unittest.TestCase):

    def test_creation(self):
        try:
            base = rnamake.base.Base()
        except:
            self.fail("got an error I did not expect")

    def test_option(self):
        base = rnamake.base.Base()
        base.option("test", 1)

        if base.option("test") != 1:
            self.fail("did not get option back")

    def test_set_options(self):
        base = rnamake.base.Base()
        options = {'test' : 10, 'test2' : "test2" }
        base.set_options(options)

        if base.option('test') != 10:
            self.fail("did not set option value correctly")


def main():
    unittest.main()

if __name__ == '__main__':
    main()
