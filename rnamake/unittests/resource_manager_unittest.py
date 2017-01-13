import unittest
from rnamake import resource_manager
import rnamake.settings as settings

from rnamake import exceptions

class ResourceManagerUnittest(unittest.TestCase):

    def setUp(self):
        self.rm = resource_manager.ResourceManager()

    def test_creation(self):
        m = self.rm.get_motif(name="HELIX.IDEAL")

    def test_get_motif(self):
        m = self.rm.get_motif(name="HELIX.IDEAL")

        with self.assertRaises(exceptions.ResourceManagerException):
            self.rm.get_motif(name="test")

    def _test_get_motif_step(self):
        m = self.rm.get_motif(end_id="GG_LL_CC_RR")

    def test_add_motif(self):
        path = settings.UNITTEST_PATH + "/resources/motifs/tetraloop_receptor_min"
        self.rm.add_motif_from_file(path)
        m = self.rm.get_motif(name="tetraloop_receptor_min",
                                 end_name="A228-A246")

        if m.name != 'tetraloop_receptor_min':
            self.fail('did not get correct motif back')



def main():
    unittest.main()

if __name__ == '__main__':
    main()
