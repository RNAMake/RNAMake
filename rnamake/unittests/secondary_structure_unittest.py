import unittest
import build
import rnamake.sqlite_library as sqlite_library
import rnamake.secondary_structure as secondary_structure
import rnamake.resource_manager as rm
import rnamake.motif_tree as motif_tree
import rnamake.motif_factory as motif_factory
import rnamake.ss_tree as ss_tree
import rnamake.secondary_structure_factory as ssfactory
import rnamake.motif_tree as motif_tree
import rnamake.setup.motif_library as motif_library
import rnamake.motif_type as motif_type
import rnamake.settings as settings
from rnamake import secondary_structure_factory as ssf

class StructureUnittest(unittest.TestCase):

    def test_creation(self):
        ss = secondary_structure.Structure(sequence="AGCU+AGCU",
                                           dot_bracket="((((+))))")

    def test_find_residue(self):
        ss = secondary_structure.Structure(sequence="AGCU+AGCU",
                                           dot_bracket="((((+))))")
        r = ss.get_residue(1, "A")
        if r is None:
            self.fail("did not find a known residue, find_residue not working")

        r2 = ss.get_residue(uuid=r.uuid)
        if r2 is None:
            self.fail("did not find a known residue, find_residue not working")

        try:
            ss.get_residue(chain_id="A")
        except ValueError:
            pass
        except:
            self.fail("did not get expected error")

    def test_copy(self):
        ss = secondary_structure.Structure(sequence="AGCU+AGCU",
                                           dot_bracket="((((+))))")
        c_ss = ss.copy()
        for r in ss.residues():
            if c_ss.get_residue(uuid=r.uuid) is None:
                self.fail("cannot find residue in copy")

    def test_to_str(self):
        ss = secondary_structure.Structure(sequence="AGCU+AGCU",
                                           dot_bracket="((((+))))")

        s = ss.to_str()
        c_ss = secondary_structure.str_to_structure(s)
        for r in ss.residues():
            if c_ss.get_residue(r.num, r.chain_id) is None:
                self.fail("cannot find residue in to_str")

    def _test_add_motif(self):
        seq1 = 'CTAGGATATGGAAGATCCTCGGGAACGAGGATCTTCCTAAGTCCTAG'
        seq2 = 'CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG'

        db1  = '(((((((..((((((((((((....))))))))))))...)))))))'
        db2  = '((((((....((((((((((((....))))))))))))....))))))'

        ss = ssfactory.factory.get_structure(seq1 +"+" + seq2, db1 + "+" + db2, to_RNA=1)

        rm.manager.add_motif('resources/motifs/GAAA_tetraloop.pdb')
        rm.manager.add_motif('resources/motifs/GGAA_tetraloop.pdb')
        m1 = rm.manager.get_motif(name='GAAA_tetraloop')
        m2 = rm.manager.get_motif(name='GGAA_tetraloop')
        ss.add_motif(m1.secondary_structure, m1.name)
        ss.add_motif(m2.secondary_structure, m2.name)
        ss_m = ss.motif('GGAA_tetraloop')
        last_end = None
        for i, end in enumerate(ss_m.ends):
            if ss_m.end_ids[i] == 'GGGAAC_LUUUUR_CCUGUGUC_LLULUULL_GAAUCUGG_RRUURURR':
                last_end = end
                break
        conn = ss.motif_topology_from_end(ss.ends[1], last_end=last_end)
        mtt = motif_tree_topology.MotifTreeTopology(conn)
        mt = motif_tree.motif_tree_from_topology(mtt, sterics=0)
        #mt.write_pdbs()


class MotifUnittest(unittest.TestCase):

    def test_creation(self):
        m = secondary_structure.Motif()
        m1 = ssf.factory.motif("AGCU+AGCU", "((((+))))")
        if len(m1.residues()) != 8:
            self.fail("did not get the correct number of residues")

        if len(m1.chains()) != 2:
            self.fail("did not get the correct number of chains")

    def test_copy(self):
        m = ssf.factory.motif("AGCU+AGCU", "((((+))))")
        m_copy = m.copy()

        for r in m.residues():
            if m_copy.get_residue(uuid=r.uuid) is None:
                self.fail("did not find residue in motif copy")

    def test_to_str(self):
        m = ssf.factory.motif("AGCU+AGCU", "((((+))))")
        #print m.to_str()
        m1 = secondary_structure.str_to_motif(m.to_str())

    def test_copy_w_res(self):
        m = ssf.factory.motif("AGCU+AGCU", "((((+))))")
        res = {r.uuid  : r for r in m.residues()}
        bps = {bp.uuid : bp for bp in m.basepairs}
        m_copy = m.copy_w_res(res, bps)

        for r in m.residues():
            if m_copy.get_residue(uuid=r.uuid) != r:
                self.fail("did not copy_w_res correctly")


class PoseUnittest(unittest.TestCase):

    def test_creation(self):
        p = secondary_structure.Pose()

        p1 = ssf.factory.pose("AGCU+AGCU", "((((+))))")
        if len(p1.motifs) != 3:
            self.fail("did not get the right number of motifs")

    def test_copy(self):
        p = ssf.factory.pose("AGCU+AGCU", "((((+))))")
        p_copy = p.copy()
        if len(p_copy.motifs) != 3:
            self.fail("did not get the right number of motifs")

    def test_to_str(self):
        p = ssf.factory.pose("AGCU+AGCU", "((((+))))")
        s = p.to_str()
        p1 = secondary_structure.str_to_pose(s)

        if len(p1.motifs) != 3:
            self.fail("did not convert from string properly")

def main():
    unittest.main()

if __name__ == '__main__':
    main()