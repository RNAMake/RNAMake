import unittest
import os
from rnamake import motif, motif_type, residue_type, settings, structure, x3dna
from rnamake import rna_structure, util, motif_state
from rnamake.primitives.rna_structure import ends_from_basepairs, assign_end_id, end_id_to_seq_and_db

import numerical, is_equal

def motif_from_pdb(pdb_path, rts):
    s = structure.structure_from_pdb(pdb_path, rts)
    bps = rna_structure.basepairs_from_x3dna(pdb_path, s)
    ends = ends_from_basepairs(s, bps)
    for end in ends:
        bps.remove(end)

    end_ids = []
    for end in ends:
        end_id = assign_end_id(s, bps, ends, end)
        end_ids.append(end_id)

    name = util.filename(pdb_path)[:-4]
    score = 0
    seq, dot_bracket = end_id_to_seq_and_db(end_ids[0])
    m = motif.Motif(s, bps, ends, end_ids, name, motif_type.HELIX, 0, dot_bracket)
    return m


class MotifUnittest(unittest.TestCase):

    def setUp(self):
        self.rts = residue_type.ResidueTypeSet()
        path = settings.RESOURCES_PATH + "base_helix/base_helix.pdb"
        self.m = motif_from_pdb(path, self.rts)

    def test_creation(self):
        self.failUnless(self.m.dot_bracket == "(((&)))")

    def test_to_str(self):
        m = self.m
        s = m.to_str()
        m2 = motif.Motif.from_str(s, self.rts)

        self.failUnless(m.name == m2.name)

    def test_copy(self):
        m = self.m
        m2 = motif.Motif.copy(m)

        self.failUnless(m.name == m2.name)

    def test_get_state(self):
        m = self.m
        ms = m.get_state()

        self.failUnless(ms.num_ends() == m.num_ends())
        self.failUnless(numerical.are_points_equal(ms.get_end(0).d,
                                                   m.get_end(0).d))

        self.failUnless(ms.num_res() == m.num_res())
        self.failUnless(ms.dot_bracket == m.dot_bracket)

        for r in m:
            self.failUnless(ms.get_residue(uuid=r.uuid) is not None)

        s = ms.to_str()
        ms_copy = motif_state.Motif.from_str(s)

        self.failUnless(ms_copy.num_res() == m.num_res())
        self.failUnless(ms_copy.mtype == m.mtype)

        ms_copy = motif_state.Motif.copy(m)

        self.failUnless(ms_copy.num_res() == m.num_res())
        self.failUnless(ms_copy.mtype == m.mtype)

    def test_state_align(self):
        m1 = self.m
        m2 = motif.Motif.copy(m1)
        ms1 = m1.get_state()
        ms2 = m2.get_state()

        motif.align_motif(m1.get_end(1), m2.get_end(0), m2)
        motif_state.align_motif_state(ms1.get_end(1), ms2)

        self.failUnless(numerical.are_points_equal(m2.get_end(1).d, ms2.get_end(1).d))
        self.failUnless(numerical.are_matrices_equal(m2.get_end(1).r, ms2.get_end(1).r))

        m3 = motif.Motif.copy(m1)
        motif.align_motif(ms2.get_end(0), m3.get_end(0), m3)

        for i in range(m3.num_res()):
            self.failUnless(is_equal.are_residues_equal(m3.get_residue(index=i),
                                                        m2.get_residue(index=i)))



def main():
    unittest.main()

if __name__ == '__main__':
    main()
