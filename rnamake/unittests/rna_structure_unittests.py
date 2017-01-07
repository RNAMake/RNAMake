import unittest
from rnamake import settings, rna_structure, residue_type

import util, instances
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

def main():
    unittest.main()

if __name__ == '__main__':
    main()
