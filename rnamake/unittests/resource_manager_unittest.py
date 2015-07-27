import unittest
import rnamake.resource_manager as rm
import rnamake.settings as settings
import rnamake.motif_tree as motif_tree

class ResourceManagerUnittest(unittest.TestCase):

    def test_creation(self):
        m = rm.manager.get_motif(name="HELIX.IDEAL")

    def test_get_motif(self):
        m = rm.manager.get_motif(name="HELIX.IDEAL")

        try:
            m = rm.manager.get_motif(name="test")
            self.fail()
        except ValueError:
            pass
        except:
            self.fail("unexpected error")

    def test_get_motif_step(self):
        m = rm.manager.get_motif(end_id="GG_LL_CC_RR")

    def test_add_motif(self):
        path = settings.UNITTEST_PATH + "/resources/motifs/tetraloop_receptor_min"
        rm.manager.add_motif(path)
        mt = motif_tree.MotifTree()
        m = rm.manager.get_motif(name="tetraloop_receptor_min",
                                 end_name="A228-A246")
        if m.name != 'tetraloop_receptor_min':
            self.fail('did not get correct motif back')


def main():
    unittest.main()

if __name__ == '__main__':
    main()
