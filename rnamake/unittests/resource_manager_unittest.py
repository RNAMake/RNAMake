import unittest
import numpy as np
from rnamake import resource_manager
import rnamake.settings as settings

from rnamake import exceptions
import numerical

class ResourceManagerUnittest(unittest.TestCase):

    def setUp(self):
        self.rm = resource_manager.ResourceManager()

    def test_get_motif(self):
        m = self.rm.get_motif(name="HELIX.IDEAL")
        self.failUnless(m.name == "HELIX.IDEAL")

        with self.assertRaises(exceptions.ResourceManagerException):
            self.rm.get_motif(name="test")

    def _test_get_motif_step(self):
        m = self.rm.get_motif(end_id="GG_LL_CC_RR")
        self.failUnless(m.get_end(0) == "GG_LL_CC_RR")

    def test_add_motif(self):
        path = settings.UNITTEST_PATH + "/resources/motifs/tetraloop_receptor_min"
        self.rm.add_motif_from_file(path)
        m = self.rm.get_motif(name="tetraloop_receptor_min",
                                 end_name="A228-A246")

        if m.name != 'tetraloop_receptor_min':
            self.fail('did not get correct motif back')

    def test_get_motif_with_new_alignment(self):
        rm = self.rm
        m = rm.get_motif(name="BP.1")
        m.move(np.array([10, 10, 10]))

        m2 = rm.get_motif_with_new_alignment(m, 1)
        self.failUnless(numerical.are_points_equal(m.get_end(0).d, m2.get_end(1).d))


def main():
    unittest.main()

if __name__ == '__main__':
    main()
