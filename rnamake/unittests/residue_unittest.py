import sys
import unittest
import logging
import numpy as np
import rnamake.residue_type
import rnamake.residue
import util


class ResidueUnittest(unittest.TestCase):

    def test_creation(self):
        gtype = rnamake.residue_type.get_rtype("GUA")
        try:
            res = rnamake.residue.Residue(gtype, "GUA", 1, "A")
        except:
            self.fail("cannot creat Residue object sucessfully")

    def test_setup_atoms_log1(self):
        gtype = rnamake.residue_type.get_rtype("GUA")
        res = rnamake.residue.Residue(gtype, "GUA", 1, "A")
        atom_names = gtype.atom_map.keys()
        atoms = [
            rnamake.atom.Atom(name, np.array([0, 1, 2]))
            for name in atom_names]
        atom1 = rnamake.atom.Atom("H1", np.array([0, 1, 2]))
        atoms.append(atom1)
        output = util.get_log_output(res.setup_atoms, atoms)
        expected = "H1 not included in <Residue('GUA1 chain A')>"
        if output != expected:
            self.fail("Did not get expected logging information")

    def test_setup_atoms_log2(self):
        gtype = rnamake.residue_type.get_rtype("GUA")
        res = rnamake.residue.Residue(gtype, "GUA", 1, "A")
        atom_names = gtype.atom_map.keys()
        atoms = [
            rnamake.atom.Atom(name, np.array([0, 1, 2]))
            for name in atom_names]
        atoms.pop()
        output = util.get_log_output(res.setup_atoms, atoms)
        print output


def basic_test_suite():
    test_str = 'test_creation test_setup_atoms_log1 test_setup_atoms_log2'
    tests = test_str.split()
    return unittest.TestSuite(map(ResidueUnittest, tests))


def main():
    if len(sys.argv) == 1:
        unittest.main()
    else:
        test_type = sys.argv[1]
        if int(test_type) == util.UnittestType.BASIC:
            unittest.TextTestRunner().run(basic_test_suite())


if __name__ == '__main__':
    main()
