import unittest
import rnamake.basepair
import rnamake.structure
import util
import numpy as np

class BasepairUnittest(unittest.TestCase):

    def setUp(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6.pdb"
        struct = util.supress_log_output(rnamake.structure.Structure,
                                         path)
        r = np.eye(3)
        bp = rnamake.basepair.Basepair(struct.get_residue(num=103),
                                       struct.get_residue(num=104),
                                       r)

        self.basepair = bp

    def test_creation(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6.pdb"
        struct = util.supress_log_output(rnamake.structure.Structure,
                                         path)
        r = np.eye(3)
        try:
            bp = rnamake.basepair.Basepair(struct.get_residue(num=103),
                                           struct.get_residue(num=104),
                                           r)
        except:
            self.fail("was not expecting an error upon creation")

    def test_residues(self):
        residues = self.basepair.residues()
        if len(residues) != 2:
            self.fail()

        if residues[0] != self.basepair.res1:
            self.fail()

    def test_partner(self):
        bp = self.basepair
        partner = bp.partner(bp.res1)
        if partner != bp.res2:
            self.fail()

        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6.pdb"
        struct = util.supress_log_output(rnamake.structure.Structure,
                                         path)
        residues = struct.residues()

        try:
            partner = bp.partner(residues[-1])
            self.fail()
        except ValueError:
            pass
        except:
            self.fail("got an error I did not expect")

    def test_state_copy(self):
        bpstate = rnamake.basepair.BasepairState(np.eye(3),
                                                 np.array([1,0,0]),
                                                 [[1,0,0],[0,1,0]])

        cbpstate = bpstate.copy()
        cbpstate.r[0][0] += 1
        if bpstate.r[0][0] == cbpstate.r[0][0]:
            self.fail()

        cbpstate.sugars[0][0] += 1
        if bpstate.sugars[0][0] == cbpstate.sugars[0][0]:
            self.fail("sugars did not deep copy")


def main():
    unittest.main()

if __name__ == '__main__':
    main()
