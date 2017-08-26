import unittest
import uuid

from rnamake import secondary_structure, exceptions
import is_equal_ss

class ResidueUnittest(unittest.TestCase):

    def test_creation(self):
        secondary_structure.Residue("G", "(", 10, "A", " ", uuid.uuid1())

    def test_copy(self):
        r =  secondary_structure.Residue("G", "(", 10, "A", " ", uuid.uuid1())
        r_copy = secondary_structure.Residue.copy(r)

        self.failUnless(is_equal_ss.are_residues_equal(r, r_copy))

        r_copy.set_name("A")
        self.failUnless(r.get_name() != r_copy.get_name())

    def test_get_str(self):
        r =  secondary_structure.Residue("G", "(", 10, "A", " ", uuid.uuid1())
        s = r.get_str()
        r_copy = secondary_structure.Residue.from_str(s)

        self.failUnless(is_equal_ss.are_residues_equal(r, r_copy, check_uuid=0))


class ChainUnittest(unittest.TestCase):

    def setUp(self):
        residues = []
        for i in range(1, 10):
            r =  secondary_structure.Residue("G", "(", i, "A", " ", uuid.uuid1())
            residues.append(r)
        self.chain = secondary_structure.Chain(residues)

    def test_copy(self):
        c = self.chain
        c_copy = secondary_structure.Chain.copy(c)

        self.failIf(len(c) != len(c_copy), "copied chain is not the right size")
        for i, r in enumerate(c):
            self.failIf(r != c_copy.get_residue(i),
                        "did not get correct uuid")

    def test_first_and_last(self):
        """test first and last functions, should return the first and last
           residues of the chain respectively.

           Also tests whether an empty chain will return exception when calling
           either
        """
        c = self.chain
        if c.get_residue(0)  != c.get_first() or \
           c.get_residue(-1) != c.get_last():
            self.fail()

        chain_2 = secondary_structure.Chain([])

        with self.assertRaises(exceptions.ChainException):
            chain_2.get_first()

        with self.assertRaises(exceptions.ChainException):
            chain_2.get_last()

    def test_to_str(self):
        c = self.chain
        s = c.get_str()
        c_copy = secondary_structure.Chain.from_str(s)


class StructureUnittest(unittest.TestCase):

    def setUp(self):
        seq = "AGCU"
        chain_1_residues = []
        chain_2_residues = []
        for i, e in enumerate(seq):
            chain_1_residues.append(
                secondary_structure.Residue(e, "(", i+1, "A", " ", uuid.uuid1()))
        for i, e in enumerate(seq):
            chain_2_residues.append(
                secondary_structure.Residue(e, ")", i+len(seq)+1, "B", " ", uuid.uuid1()))

        c1 = secondary_structure.Chain(chain_1_residues)
        c2 = secondary_structure.Chain(chain_2_residues)
        all_res, chain_cuts = secondary_structure.chains_to_res_and_cuts([c1, c2])
        self.s = secondary_structure.Structure(all_res, chain_cuts)

    def test_creation(self):
        pass

    def test_find_residue(self):
        ss = self.s
        r = ss.get_residue(1, "A")
        if r is None:
            self.fail("did not find a known residue, find_residue not working")

        r2 = ss.get_residue(uuid=r.get_uuid())
        if r2 is None:
            self.fail("did not find a known residue, find_residue not working")

        with self.assertRaises(exceptions.StructureException):
            ss.get_residue(chain_id="A")

    def test_copy(self):
        ss = self.s
        c_ss = secondary_structure.Structure.copy(ss)
        for r in ss:
            if c_ss.get_residue(uuid=r.get_uuid()) is None:
                self.fail("cannot find residue in copy")

    def test_get_str(self):
        ss = self.s
        s = ss.get_str()
        c_ss = secondary_structure.Structure.from_str(s)
        for r in ss:
            if c_ss.get_residue(r.get_num(), r.get_chain_id()) is None:
                self.fail("cannot find residue in to_str")

    def test_seq_and_ss(self):
        self.failUnless(self.s.get_sequence() == "AGCU&AGCU")


class BasepairUnittest(unittest.TestCase):

    def setUp(self):
        pass


class RNAStructureUnittest(unittest.TestCase):
    def setUp(self):
        seq = "AGCU"
        chain_1_residues = []
        chain_2_residues = []
        for i, e in enumerate(seq):
            chain_1_residues.append(secondary_structure.Residue(e, "(", i+1, "A"))
        for i, e in enumerate(seq):
            chain_2_residues.append(secondary_structure.Residue(e, ")", i+len(seq)+1, "B"))

        c1 = secondary_structure.Chain(chain_1_residues)
        c2 = secondary_structure.Chain(chain_2_residues)
        s = secondary_structure.Structure([c1, c2])
        bps_indexes = [
            [1, 8],
            [2, 7],
            [3, 6],
            [4, 5]]

        bps = []
        for r_i1, r_i2 in bps_indexes:
            r1 = s.get_residue(r_i1)
            r2 = s.get_residue(r_i2)
            bp = secondary_structure.Basepair(r1.uuid, r2.uuid, basepair.calc_bp_name([r1, r2]))
            bps.append(bp)

        ends = rna_structure.ends_from_basepairs(s, bps)
        end_id_1 = rna_structure.assign_end_id(s, bps, ends[0])
        end_id_2 = rna_structure.assign_end_id(s, bps, ends[1])
        end_ids = [end_id_1, end_id_2]
        self.rna_struc = secondary_structure.RNAStructure(s, bps, ends, end_ids)

    def test_creation(self):
        pass

    def test_copy(self):
        rna_struc_copy = secondary_structure.RNAStructure.copy(self.rna_struc)
        self.failUnless(self.rna_struc.get_end(0) == rna_struc_copy.get_end(0))

    def test_to_str(self):
        pass


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

    def test_get_helices(self):
        p = ssf.factory.pose("AGCUAGG+CCAGCU",
                             "((((.((+))))))")
        p.build_helices()
        self.failIf(len(p.helices) != 2, "did not build the right number of helices")

def main():
    unittest.main()

if __name__ == '__main__':
    main()