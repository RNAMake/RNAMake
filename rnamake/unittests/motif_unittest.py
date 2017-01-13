import unittest
import os
from rnamake import motif, motif_type, residue_type, settings, structure, x3dna
from rnamake import rna_structure, util, motif_state
from rnamake.primitives.rna_structure import ends_from_basepairs, assign_end_id, end_id_to_seq_and_db

import numerical

def motif_from_pdb(pdb_path, rts):
    s = structure.structure_from_pdb(pdb_path, rts)
    bps = rna_structure.basepairs_from_x3dna(pdb_path, s)
    ends = ends_from_basepairs(s, bps)
    end_ids = []
    for end in ends:
        end_id = assign_end_id(s, bps, end)
        end_ids.append(end_id)
    name = util.filename(pdb_path)[:-4]
    score = 0
    seq, dot_bracket = end_id_to_seq_and_db(end_ids[0])
    m = motif.Motif(s, bps, ends, end_ids, name, motif_type.HELIX, 0,
                    dot_bracket=dot_bracket)
    return m


class MotifUnittest(unittest.TestCase):

    def setUp(self):
        self.rts = residue_type.ResidueTypeSet()
        path = settings.RESOURCES_PATH + "base_helix/base_helix.pdb"
        self.m = motif_from_pdb(path, self.rts)

    def test_creation(self):
        pass

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

        for r in m.iter_res():
            self.failUnless(ms.get_residue(uuid=r.uuid) is not None)

        s = ms.to_str()
        ms_copy = motif_state.Motif.from_str(s)

        self.failUnless(ms_copy.num_res() == m.num_res())
        self.failUnless(ms_copy.mtype == m.mtype)

        ms_copy = motif_state.Motif.copy(m)

        self.failUnless(ms_copy.num_res() == m.num_res())
        self.failUnless(ms_copy.mtype == m.mtype)


    def _test_state(self):
        ms1 = rm.manager.get_state(name="HELIX.IDEAL.2")
        ms2 = rm.manager.get_state(name="HELIX.IDEAL.2")

        motif.align_motif_state(ms1.end_states[1], ms2)

        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")

        motif.align_motif(m1.ends[1].state(), m2.ends[0], m2)

        if util.distance(ms2.end_states[1].d, m2.ends[1].d()) > 0.01:
            print ms2.end_states[1].d
            print m2.ends[1].d()
            self.fail("motif state did not act like a motif for origin")

        if util.matrix_distance(ms2.end_states[1].r, m2.ends[1].r()) > 0.01:
            print ms2.end_states[1].r
            print m2.ends[1].r()
            self.fail("motif state did not act like a motif for rotation")

    def _test_to_str(self):
        m = self.motif
        s = m.to_str()
        m1 = rnamake.motif.str_to_motif(s)
        if len(m1.residues()) != 10:
            self.fail("did not copy all residues correctly")

    """def test_get_end_id(self):
        m = rm.manager.get_motif(name="HELIX.IDEAL")
        end_id = m.end_index_with_id('GG_LL_CC_RR')

    def _test_protein_beads(self):
        path = files.GROUP_2_INTRON_PDB_PATH
        m = rnamake.motif_factory.factory.motif_from_file(path, include_protein=1)

        beads = m.protein_beads
        #basic_io.beads_to_pdb("test.pdb", beads)"""


def main():
    unittest.main()

if __name__ == '__main__':
    main()
