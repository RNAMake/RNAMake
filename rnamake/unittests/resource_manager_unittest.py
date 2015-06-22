import unittest
import rnamake.resource_manager

class ResourceManagerUnittest(unittest.TestCase):

    def test_creation(self):
        m = rnamake.resource_manager.manager.get_motif("HELIX.IDEAL")

    def test_get_motif(self):
        rm = rnamake.resource_manager.ResourceManager()
        m = rm.get_motif("HELIX.IDEAL")

        try:
            m = rm.get_motif("test")
            self.fail()
        except ValueError:
            pass
        except:
            self.fail("unexpected error")

    def test_get_motif_step(self):
        rm = rnamake.resource_manager.ResourceManager()
        m = rm.get_motif("GC=GC")



def main():
    unittest.main()

if __name__ == '__main__':
    main()
