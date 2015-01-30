import unittest
import rnamake.basepair
import rnamake.structure
import util
import numpy as np
import copy


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
                                                 np.array([1, 0, 0]),
                                                 [[1, 0, 0], [0, 1, 0]])

        cbpstate = bpstate.copy()
        cbpstate.r[0][0] += 1
        if bpstate.r[0][0] == cbpstate.r[0][0]:
            self.fail()

        cbpstate.sugars[0][0] += 1
        if bpstate.sugars[0][0] == cbpstate.sugars[0][0]:
            self.fail("sugars did not deep copy")

    def test_state(self):
        bp = self.basepair
        old_state = bp.bp_state
        old_d = old_state.d
        old_sugars = copy.deepcopy(old_state.sugars)
        bp.res1.get_atom("C1'").coords[0] += 100
        new_d = old_state.d
        if old_d[0] != new_d[0]:
            self.fail("state changed")

        new_state = bp.state()
        if new_state.d[0] == old_d[0]:
            self.fail("did not account to change coords")

    def test_copy(self):
        bp = self.basepair
        cbp = bp.copy()
        old_state = cbp.bp_state
        old_sugars = copy.deepcopy(old_state.sugars)
        cbp.res1.get_atom("C1'").coords[0] += 100

        if old_sugars[0][0] == old_state.sugars[0][0]:
            self.fail("failed to update sugars after copy")

    def test_pdb_str(self):
        bp = self.basepair
        s = bp.to_pdb_str()


def main():
    unittest.main()

if __name__ == '__main__':
    main()
