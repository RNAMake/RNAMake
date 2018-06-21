import unittest
import os
import rnamake.resource_manager as rm
import rnamake.motif as motif
import rnamake.settings
import rnamake.motif_type
import rnamake.transform
import rnamake.motif_factory
import rnamake.util as util
import numerical
import numpy as np
from rnamake import secondary_structure_factory as ssf
from rnamake import transformations
from rnamake import basic_io

import files, instances


class MotifUnittest(unittest.TestCase):

    def setUp(self):
        self.motif = instances.motif()

    def test_state_1(self):
        ms1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        state = ms1.get_state()

        s = state.to_str()
        state2 = motif.str_to_motif_state(s)

    def test_state(self):
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

    def test_get_basepair(self):
        m = self.motif
        bp = m.basepairs[0]

        found = m.get_basepair(res1=bp.res1, res2=bp.res2)
        if bp != found[0]:
            self.fail("did not retreive correct basepair")

        found = m.get_basepair(uuid1=bp.res1.uuid, uuid2=bp.res2.uuid)
        if bp != found[0]:
            self.fail("did not retreive correct basepair")

    def test_get_beads(self):
        m = self.motif
        beads = m.get_beads()
        org_count = len(beads)

        beads = m.get_beads(excluded_res=[m.structure.residues()[1]])
        diff = org_count - len(beads)
        if diff != 3:
            self.fail("did not exclude res properly")

        beads = m.get_beads([m.basepairs[0]])
        diff = org_count - len(beads)
        if diff != 6:
            self.fail("did not exclude ends properly")

    def test_to_str(self):
        m = self.motif
        s = m.to_str()
        m1 = rnamake.motif.str_to_motif(s)
        if len(m1.residues()) != 10:
            self.fail("did not copy all residues correctly")

    def test_copy(self):
        m = self.motif
        cm = m.copy()

    def test_transform(self):
        m = self.motif
        old_r = m.basepairs[0].state().r
        r = transformations.random_rotation_matrix()[:3,:3]
        d = np.random.random([3])
        t = rnamake.transform.Transform(r, d)
        m.transform(t)
        new_r = m.basepairs[0].state().r
        if numerical.are_matrices_equal(old_r, new_r):
            self.fail("rotations should be different")

    def test_get_end_id(self):
        m = rm.manager.get_motif(name="HELIX.IDEAL")
        end_id = m.end_index_with_id('GG_LL_CC_RR')

    def _test_protein_beads(self):
        path = files.GROUP_2_INTRON_PDB_PATH
        m = rnamake.motif_factory.factory.motif_from_file(path, include_protein=1)

        beads = m.protein_beads
        #basic_io.beads_to_pdb("test.pdb", beads)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
