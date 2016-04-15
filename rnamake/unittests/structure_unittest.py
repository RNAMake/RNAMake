import unittest
import rnamake.structure
import rnamake.transform
import rnamake.motif_factory
import rnamake.io

from rnamake import exceptions

import is_equal
import numpy as np

class StructureUnittest(unittest.TestCase):

    def setUp(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6.pdb"
        self.structure = rnamake.structure.structure_from_pdb(path)

    def test_creation(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6.pdb"
        rnamake.structure.structure_from_pdb(path)

    def test_build_chains_all(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/build_chains.dat"

        f = open(path)
        lines = f.readlines()
        f.close()

        skip_structs = "3SLQ".split()

        path = "/Users/josephyesselman/projects/REDESIGN/resources/non-redundant-rnas"
        for l in lines:
            spl = l.split(",")
            #print spl[0]
            if spl[0] in skip_structs:
                continue
            pdb_path = path + "/" + spl[0] + "/" + spl[0] + ".pdb"
            struct = rnamake.structure.structure_from_pdb(pdb_path)
            for ckey in spl[1:-1]:
                found = 0
                for c in struct.chains:
                    start = c.first()
                    end = c.last()
                    start_key = start.name + " " + str(start.num) + " " + start.chain_id
                    end_key = end.name + " " + str(end.num) + " " + end.chain_id
                    key = start_key + " " + end_key
                    if key == ckey:
                        found = 1
                        break
                if not found:
                    self.fail()

    def test_get_residue(self):
        struct = self.structure

        with self.assertRaises(exceptions.StructureException):
            struct.get_residue()

        res = struct.get_residue(num=107)
        if res is None:
            self.fail("should not of gotten an error")

        res = struct.get_residue(num=1000)
        if res is not None:
            self.fail("should not of gotten an error")

    def test_residues(self):
        struct = self.structure
        residues = struct.residues()
        if len(residues) != 157:
            self.fail()

    def test_atoms(self):
        struct = self.structure
        atoms = struct.atoms()
        if len(atoms) != 3357:
            self.fail()

    def test_to_str(self):
        struct = self.structure

        s = struct.to_str()
        struct_new = rnamake.io.str_to_structure(s)
        if len(struct_new.residues()) != len(struct.residues()):
            self.fail("did not get back all residues")

        if not is_equal.are_structure_equal(struct, struct_new, check_uuid=0):
            self.fail("did not produce the same structure from str")

    def test_get_beads(self):
        struct = self.structure
        beads = struct.get_beads()
        if len(beads) != 470:
            self.fail("got wrong number of beads")

        r = struct.get_residue(num=106)
        beads = struct.get_beads(excluded_res=[r])
        if len(beads) != 467:
            self.fail("got wrong number of beads")

    def test_transform(self):
        r = np.random.random([3,3])
        d = np.random.random([3])
        t = rnamake.transform.Transform(r, d)
        s = self.structure.copy()
        struct = self.structure
        s.transform(t)
        if is_equal.are_structure_equal(s, struct):
            self.fail("did not transform")

        struct.transform(t)

        if not is_equal.are_structure_equal(s, struct):
            self.fail("structures should be the same now")

    def test_move(self):
        path = "/Users/josephyesselman/projects/REDESIGN/redesign/tests/p4p6"
        d = np.random.random([3])
        s = self.structure.copy()
        struct = self.structure
        struct.move(d)
        if is_equal.are_structure_equal(s, struct):
            self.fail("did not move")

        s.move(d)
        if not is_equal.are_structure_equal(s, struct):
            self.fail("did not move")

    def test_copy(self):
        s = self.structure
        s_copy = s.copy()

        if not is_equal.are_structure_equal(s, s_copy):
            self.fail("copying did not return an indentical structure")

def main():
    unittest.main()

if __name__ == '__main__':
    main()
