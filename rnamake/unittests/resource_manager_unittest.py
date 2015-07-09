import unittest
import rnamake.resource_manager
import rnamake.settings as settings
import rnamake.motif_tree as motif_tree

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
        m = rm.get_motif("GG_LL_CC_RR")

    def test_add_motif(self):
        path = settings.UNITTEST_PATH + "/resources/motifs/tetraloop_receptor_min"
        rnamake.resource_manager.manager.add_motif(path)
        mt = motif_tree.MotifTree()
        m = rnamake.resource_manager.manager.get_motif("tetraloop_receptor_min",
                                                       "A228-A246")
        mt.add_motif(m)




def main():
    unittest.main()

if __name__ == '__main__':
    main()
