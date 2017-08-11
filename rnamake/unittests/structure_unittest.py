import unittest
import warnings
import rnamake.transform

from rnamake import all_atom, exceptions, residue_type, motif_state
from rnamake import settings
from instances import transform_instances

import is_equal, instances
import numpy as np

#warnings.simplefilter("ignore")

class StructureUnittest(unittest.TestCase):

    def setUp(self):
        path = settings.UNITTEST_PATH + "resources/p4p6.pdb"
        self.rts = residue_type.ResidueTypeSet()
        self.structure = all_atom.structure_from_pdb(path, self.rts)

    def test_creation(self):
        path = settings.UNITTEST_PATH + "resources/p4p6.pdb"
        all_atom.structure_from_pdb(path, self.rts)

    def test_new_chains(self):
        s = self.structure
        chains = s.get_chains()
        new_chains = []
        new_res = []
        all_res = []
        chain_cuts = []
        i = 0
        for r in chains[0]:
            new_res.append(r)
            all_res.append(r)
            if i > 10:
                new_chain = all_atom.Chain(new_res)
                new_chains.append(new_chain)
                new_res = []
                chain_cuts.append(len(all_res))
                i = 0
            i += 1

        new_chains.append(all_atom.Chain(new_res))
        chain_cuts.append(len(all_res))
        new_s = all_atom.Structure(all_res, chain_cuts)
        new_s_chains = new_s.get_chains()

        for i in range(len(new_chains)):
            self.failUnless(is_equal.are_chains_equal(new_chains[i], new_s_chains[i]))

    # TODO move to integration
    def _test_build_chains_all(self):
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

    def test_iter_residues(self):
        s = self.structure
        res = []
        for r in s:
            res.append(r)
        self.failUnless(len(res) == s.get_num_residues())

    def test_get_str(self):
        struct = self.structure

        s = struct.get_str()
        struct_new = all_atom.Structure.from_str(s, self.rts)
        if struct_new.get_num_residues() != struct.get_num_residues():
            self.fail("did not get back all residues")

        if not is_equal.are_structure_equal(struct, struct_new, check_uuid=0):
            self.fail("did not produce the same structure from str")

    def test_transform(self):
        r = np.random.random([3,3])
        d = np.random.random([3])
        t = rnamake.transform.Transform(r, d)
        s = all_atom.Structure.copy(self.structure)
        struct = self.structure
        s.transform(t)
        if is_equal.are_structure_equal(s, struct):
            self.fail("did not transform")

        struct.transform(t)

        if not is_equal.are_structure_equal(s, struct):
            self.fail("structures should be the same now")

        s2 = all_atom.Structure.copy(self.structure)
        s3 = all_atom.Structure.copy(self.structure)
        s2.transform(transform_instances.transform_indentity())

        if not is_equal.are_structure_equal(s2, s3):
            self.fail("did not transform")

    def test_move(self):
        path = "/Users/josephyesselman/projects/REDESIGN/redesign/tests/p4p6"
        d = np.random.random([3])
        s = all_atom.Structure.copy(self.structure)
        struct = self.structure
        struct.move(d)
        if is_equal.are_structure_equal(s, struct):
            self.fail("did not move")

        s.move(d)
        if not is_equal.are_structure_equal(s, struct):
            self.fail("did not move")

    def test_copy(self):
        s = self.structure
        s_copy = all_atom.Structure.copy(s)

        if not is_equal.are_structure_equal(s, s_copy):
            self.fail("copying did not return an indentical structure")

    def test_get_state(self):
        s = self.structure
        ss = s.get_state()

        self.failUnless(len(s) == len(ss))

        ss_copy = motif_state.Structure.copy(ss)

        s = ss.get_str()
        ss_copy = motif_state.Structure.from_str(s)

        self.failUnless(len(self.structure) == len(ss_copy))


def main():
    unittest.main()

if __name__ == '__main__':
    main()
