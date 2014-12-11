import sys
import unittest
import logging
import numpy as np
import rnamake.residue_type
import rnamake.residue
import rnamake.io
import rnamake.settings
import rnamake.util
import util


class ResidueUnittest(unittest.TestCase):

    def setUp(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/res_strs.dat"
        f = open(path)
        lines = f.readlines()
        f.close()

        residues = []
        logging.disable(60)
        for l in lines:
            res = rnamake.io.str_to_residue(l)
            residues.append(res)
        logging.disable(0)
        self.residues = residues

    def test_creation(self):
        """
        make sure creating a Residue object does not throw any exceptions
        """
        gtype = rnamake.residue_type.get_rtype("GUA")
        try:
            res = rnamake.residue.Residue(gtype, "GUA", 1, "A")
        except:
            self.fail("cannot creat Residue object sucessfully")

    def test_setup_atoms_log1(self):
        """
        first tests in logging, make sure that
        """
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
        expected = "O3' is undefined in <Residue('GUA1 chain A')>"
        if output != expected:
            self.fail()

    def test_get_atom(self):
        gtype = rnamake.residue_type.get_rtype("GUA")
        res = rnamake.residue.Residue(gtype, "GUA", 1, "A")
        atom_names = gtype.atom_map.keys()
        atoms = [
            rnamake.atom.Atom(name, np.array([0, 1, 2]))
            for name in atom_names]
        res.setup_atoms(atoms)

        # can find a real atom
        p_atom = res.get_atom("P")
        self.assertIs(p_atom, res.atoms[0])

        # should throw and error if no error exist
        try:
            no_atom = res.get_atom("P1")
        except KeyError:
            pass
        except:
            self.fail("did not expect this error")

    def test_str_to_residue(self):
        """
        test
        """
        path = rnamake.settings.UNITTEST_PATH + "resources/res_strs.dat"
        f = open(path)
        lines = f.readlines()
        f.close()

        try:
            logging.disable(60)
            res = rnamake.io.str_to_residue(lines[0])
            logging.disable(0)
        except:
            self.fail()

    def test_connected_to(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/res_strs.dat"
        f = open(path)
        lines = f.readlines()
        f.close()

        logging.disable(60)
        res1 = rnamake.io.str_to_residue(lines[0])
        res2 = rnamake.io.str_to_residue(lines[1])
        res3 = rnamake.io.str_to_residue(lines[2])
        logging.disable(0)

        result = res1.connected_to(res2)
        if result != 1:
            self.fail()

        result = res2.connected_to(res1)
        if result != -1:
            self.fail()

        result = res1.connected_to(res3)
        if result != 0:
            self.fail()

    def test_get_beads(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/res_strs.dat"
        f = open(path)
        lines = f.readlines()
        f.close()

        residues = []
        logging.disable(60)
        for l in lines:
            res = rnamake.io.str_to_residue(l)
            residues.append(res)
        logging.disable(0)

        for r in residues:
            beads = r.get_beads()
            if r.get_atom("P") is None and len(beads) != 2:
                self.fail()
            elif r.get_atom("P") is not None and len(beads) != 3:
                self.fail()

        # check calculations
        phos_atom_names = "P OP2 OP1".split(" ")
        sugar_atom_names = "O5' C5' C4' O4' C3' O3' C1' C2' O2'".split(" ")
        phos_atoms, sugar_atoms, base_atoms = [], [], []
        for a in residues[1].atoms:
            if a.name in phos_atom_names:
                phos_atoms.append(a)
            elif a.name in sugar_atom_names:
                sugar_atoms.append(a)
            else:
                base_atoms.append(a)

        beads = residues[1].get_beads()
        phos_center = rnamake.util.center(phos_atoms)
        sugar_center = rnamake.util.center(sugar_atoms)
        base_center = rnamake.util.center(base_atoms)

        self.assertAlmostEqual(
            rnamake.util.distance(
                phos_center,
                beads[0].center),
            0)

        self.assertAlmostEqual(
            rnamake.util.distance(
                sugar_center,
                beads[1].center),
            0)

        self.assertAlmostEqual(
            rnamake.util.distance(
                base_center,
                beads[2].center),
            0)

    def test_copy(self):
        res = self.residues[0]
        copy_res = res.copy()
        copy_res.name = "test"
        if res.name == copy_res.name:
            self.fail("did not copy name correctly")
        copy_res.num = 1000
        if res.num == copy_res.num:
            self.fail("did not copy num correctly")

        print copy_res.name, res.name

    def test_new_uuid(self):
        res = self.residues[0]
        old_uuid = res.uuid
        res.new_uuid()
        if old_uuid == res.uuid:
            self.fail("did not assign new uuid to res")

    def test_to_pdb_str(self):
        res = self.residues[0]
        s = res.to_pdb_str()


def main():
    if len(sys.argv) == 1:
        unittest.main()
    else:
        test_type = sys.argv[1]
        if int(test_type) == util.UnittestType.BASIC:
            unittest.TextTestRunner().run(basic_test_suite())


if __name__ == '__main__':
    main()
