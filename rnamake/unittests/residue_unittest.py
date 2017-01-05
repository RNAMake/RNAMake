import sys
import unittest
import numpy as np

from rnamake.residue import Residue, Atom, Bead, BeadType
from rnamake import residue_type, settings, util, exceptions

import is_equal


class ResidueUnittest(unittest.TestCase):

    def setUp(self):
        self.rts = residue_type.ResidueTypeSet()

        path = settings.UNITTEST_PATH + "resources/res_strs.dat"
        f = open(path)
        lines = f.readlines()
        f.close()

        residues = []
        for l in lines:
            res = Residue.from_str(l, self.rts)
            residues.append(res)
        self.residues = residues

    def test_creation(self):
        """
        make sure creating a Residue object does not throw any exceptions
        """

        try:
            #also check name changing works
            rtype = self.rts.get_type("GUA")
            atoms = [Atom("O1P", np.array([1, 2, 3]))]
            res = Residue(atoms, rtype, "GUA", 1, "A")
        except:
            self.fail("cannot creat Residue object sucessfully")

    def test_get_atom(self):
        gtype = self.rts.get_type("GUA")
        atoms = [
            Atom("P", np.array([0, 1, 2]))
        ]
        res = Residue(atoms, gtype, "GUA", 1, "A")

        # can find a real atom
        p_atom = res.get_atom("P")
        self.assertIs(p_atom, res.get_atom(index=0))

        # should throw and error if no error exist
        with self.assertRaises(exceptions.ResidueException):
            res.get_atom("P1")

        # invalid num
        with self.assertRaises(exceptions.ResidueException):
            res.get_atom(index=999)

        # valid name but not initialized
        with self.assertRaises(exceptions.ResidueException):
            res.get_atom("OP1")

        with self.assertRaises(exceptions.ResidueException):
            res.get_atom(index=1)

    def test_has_atom(self):
        gtype = self.rts.get_type("GUA")
        atoms = [
            Atom("P", np.array([0, 1, 2]))
        ]
        res = Residue(atoms, gtype, "GUA", 1, "A")

        self.failUnless(res.has_atom("P") == True)
        self.failUnless(res.has_atom("P1") == False)
        self.failUnless(res.has_atom("OP1") == False )
        self.failUnless(res.has_atom(index=0) == True)
        self.failUnless(res.has_atom(index=999) == False)

        res = self.residues[0]
        # first residue has no P
        self.failUnless(res.has_atom("P") == False)
        self.failUnless(res.has_atom("C1'") == True)

    def test_to_str(self):
        res = self.residues[0]
        s = res.to_str()
        res2 = Residue.from_str(s, self.rts)

        self.failUnless(is_equal.are_residues_equal(res, res2, check_uuid=0))

    def test_connected_to(self):
        res1 = self.residues[0]
        res2 = self.residues[1]
        res3 = self.residues[2]

        self.failUnless(res1.connected_to(res2) == 1)
        self.failUnless(res2.connected_to(res1) == -1)
        self.failUnless(res1.connected_to(res3) == 0)

    def test_beads(self):
        res1 = self.residues[0]
        res2 = self.residues[10]

        for b1 in res1.iter_beads():
            for b2 in res2.iter_beads():
                self.failIf(b1.distance(b2) < settings.CLASH_RADIUS)

    def test_copy(self):
        res = self.residues[0]
        copy_res = Residue.copy(res)

        self.failUnless(is_equal.are_residues_equal(res, copy_res, check_uuid=0))

    def test_move(self):
        r = self.residues[0]
        c1 = r.center()

        r.move(np.array([1, 0, 0]))
        c2 = r.center()

        self.failUnlessAlmostEqual(1.0, util.distance(c1, c2))


class BeadUnittest(unittest.TestCase):

    def setUp(self):
        self.b = Bead(np.array([1,1,1]), BeadType.PHOS)

    def test_creation(self):

        with self.assertRaises(exceptions.ResidueException):
            Bead(np.array([1,1,1]), 4)

    def test_copy(self):
        b_copy = Bead.copy(self.b)
        d = b_copy.distance(self.b)
        self.failUnless(d < 0.0001)

    def test_type_name(self):
        name = self.b.type_name()
        if name != "PHOSPHATE":
            self.fail("did not get the right bead name")

def main():
    unittest.main()


if __name__ == '__main__':
    main()
