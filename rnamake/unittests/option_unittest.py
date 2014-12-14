import unittest
import rnamake.option

class OptionUnittest(unittest.TestCase):

    def test_creation(self):
        option = rnamake.option.Option(5, 'float')
        options = rnamake.option.Options()
        options.add("test", 5)

    def test_add(self):
        options = rnamake.option.Options()
        options.add("test", 5,)

        try:
            options.add("test", 5)
            self.fail()
        except ValueError:
            pass
        except:
            self.fail("got an error I did not expect")

    def test_set(self):
        options = rnamake.option.Options()
        options.add("test", 5)
        #should work
        options.set("test", 10)

        try:
            options.set("test2", 10)
        except ValueError:
            pass
        except:
            self.fail("did not catch this error")

        try:
            options.set("test", 10.0)
        except ValueError:
            pass
        except:
            self.fail("did not catch this error")

        try:
            options.set("test", "test")
        except ValueError:
            pass
        except:
            self.fail("did not catch this error")

    def test_add(self):
        options = rnamake.option.Options()
        options.add("test", 5)
        if options.get("test") != 5:
            self.fail("did not retreive value correctly")

    def test_contains(self):
        options = rnamake.option.Options()
        options.add("test", 5)
        result = "test" in options
        if result != True:
            self.fail()

def main():
    unittest.main()

if __name__ == '__main__':
    main()
