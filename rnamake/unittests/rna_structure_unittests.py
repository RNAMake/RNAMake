import unittest
from rnamake import settings, rna_structure, residue_type, basepair, util
from rnamake.primitives.rna_structure import align_rna_structure

import instances
import numpy as np
import copy
import random

from instances.transform_instances import transform_random
import numerical


class RNAStructureUnittest(unittest.TestCase):

    def setUp(self):
        self.rts = residue_type.ResidueTypeSet()
        path = settings.MOTIF_DIRS + "helices/HELIX.IDEAL/HELIX.IDEAL.pdb"
        self.rna_struc = rna_structure.rna_structure_from_pdb(path, self.rts)

    def test_creation(self):
        rs = self.rna_struc
        self.failUnless(rs.sequence() == "GG&CC")
        self.failUnless(rs.dot_bracket() == "((&))")

    def test_copy(self):
        rs = self.rna_struc
        rs_copy = rna_structure.RNAStructure.copy(rs)

        for r in rs.iter_res():
            self.failUnless(rs_copy.get_residue(uuid=r.uuid) is not None)

        for bp in rs.iter_basepairs():
            self.failUnless(rs_copy.get_basepair(bp_uuid=bp.uuid) is not None)

        rs_copy_2 = rna_structure.RNAStructure.copy(rs, new_uuid=1)
        for r in rs.iter_res():
            self.failUnless(rs_copy_2.get_residue(uuid=r.uuid) is None)

        for bp in rs.iter_basepairs():
            self.failUnless(rs_copy_2.get_basepair(bp_uuid=bp.uuid) is None)

    def test_to_str(self):
        rs = self.rna_struc
        s = rs.to_str()
        rs_copy = rna_structure.RNAStructure.from_str(s, self.rts)

        self.failUnless(rs.num_res() == rs_copy.num_res())
        self.failUnless(rs.num_chains() == rs_copy.num_chains())
        self.failUnless(rs.num_basepairs() == rs_copy.num_basepairs())

    def test_align(self):
        rs = self.rna_struc
        rs2 = rna_structure.RNAStructure.copy(rs)
        rs3 = rna_structure.RNAStructure.copy(rs)

        align_rna_structure(rs.get_end(1), rs2.get_end(0), rs2)
        align_rna_structure(rs2.get_end(1), rs3.get_end(0), rs3)

        # actual coords are probably aligned
        c1 = basepair.calc_center(rs.get_bp_res(rs.get_end(1)))
        c2 = basepair.calc_center(rs2.get_bp_res(rs2.get_end(0)))
        diff = util.distance(c1, c2)
        self.failUnless(diff < 0.01)

        c1 = basepair.calc_center(rs2.get_bp_res(rs2.get_end(1)))
        c2 = basepair.calc_center(rs3.get_bp_res(rs3.get_end(0)))
        diff = util.distance(c1, c2)
        self.failUnless(diff < 0.01)

        self.failUnless(util.matrix_distance(rs.get_end(1).r, rs2.get_end(0).r) < 0.1)
        self.failUnless(util.distance(rs.get_end(1).d, rs2.get_end(0).d) < 0.1)

        #rs.to_pdb("test.pdb")
        #rs2.to_pdb("test2.pdb")
        #rs3.to_pdb("test3.pdb")

    def test_get_basepair(self):
        rs = self.rna_struc
        bp = rs.get_end(0)

        found = rs.get_basepair(uuid1=bp.res1_uuid, uuid2=bp.res2_uuid)
        if bp != found:
            self.fail("did not retreive correct basepair")

    def test_beads_generation(self):
        rs = self.rna_struc
        # should not have a beads as it is the alignment end
        bp_res = rs.get_bp_res(rs.get_end(0))
        self.failUnless(bp_res[0].num_beads() != 0)
        self.failUnless(bp_res[1].num_beads() != 0)

        # an end to add too should not have beads
        bp_res = rs.get_bp_res(rs.get_end(1))
        self.failUnless(bp_res[0].num_beads() == 0)
        self.failUnless(bp_res[1].num_beads() == 0)

    def test_transform(self):
        rs = self.rna_struc
        t = transform_random()

        old_r = rs.get_end(0).r
        rs.transform(t)
        new_r = rs.get_end(0).r
        self.failUnless(numerical.are_matrices_equal(old_r, new_r) == 0)

def main():
    unittest.main()

if __name__ == '__main__':
    main()
