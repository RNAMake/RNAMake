import unittest
import uuid
import instances


from rnamake import secondary_structure, motif_tree, exceptions

from rnamake import secondary_structure_factory as ssf
from rnamake import resource_manager as rm
import rnamake.secondary_structure_factory as ssfactory

class ResidueUnittest(unittest.TestCase):

    def test_creation(self):
        secondary_structure.Residue("G", "(", 10, "A", uuid.uuid1())

    def test_copy(self):
        r =  secondary_structure.Residue("G", "(", 10, "A", uuid.uuid1())
        r_copy = r.copy()

        r.dot_bracket = ")"
        self.failIf(r.dot_bracket == r_copy.dot_bracket, "should not be equal")

    def test_to_str(self):
        r =  secondary_structure.Residue("G", "(", 10, "A", uuid.uuid1())
        s = r.to_str()
        r_copy = secondary_structure.str_to_residue(s)
        self.failIf(r.name != r_copy.name)
        self.failIf(r.dot_bracket != r_copy.dot_bracket)


class ChainUnittest(unittest.TestCase):

    def test_creation(self):
        secondary_structure.Chain()

    def test_copy(self):
        c = instances.secondary_structure_chain()
        c_copy = c.copy()

        self.failIf(len(c) != len(c_copy), "copied chain is not the right size")
        for i, r in enumerate(c):
            self.failIf(r.uuid != c_copy.residues[i].uuid,
                        "did not get correct uuid")

    def test_first_and_last(self):
        """test first and last functions, should return the first and last
           residues of the chain respectively.

           Also tests whether an empty chain will return exception when calling
           either
        """
        c = instances.secondary_structure_chain()
        if c.residues[0]  != c.first() or \
           c.residues[-1] != c.last():
            self.fail()

        chain_2 = secondary_structure.Chain()

        with self.assertRaises(exceptions.SecondaryStructureException):
            chain_2.first()

        with self.assertRaises(exceptions.SecondaryStructureException):
            chain_2.last()

    def test_to_str(self):
        r =  secondary_structure.Residue("G", "(", 10, "A", uuid.uuid1())
        s = r.to_str()
        r_copy = secondary_structure.str_to_residue(s)


class StructureUnittest(unittest.TestCase):

    def test_creation(self):
        secondary_structure.Structure(sequence="AGCU+AGCU",
                                      dot_bracket="((((+))))")

        with self.assertRaises(exceptions.SecondaryStructureException):
            secondary_structure.Structure(sequence="AGCU+AGCU",
                                          dot_bracket="^((((+))))")

        #with self.assertRaises(exceptions.SecondaryStructureException):
        #    secondary_structure.Structure(sequence="KGCU+AGCU",
        #                                  dot_bracket="((((+))))")

        with self.assertRaises(exceptions.SecondaryStructureException):
            secondary_structure.Structure(sequence="GCU+AGCU",
                                          dot_bracket="((((+))))")


        seq = ""
        db = ""
        for i in range(100):
            seq += "A&"
            db += ".A"

        #recycles chain ids does not run out of letters
        secondary_structure.Structure(sequence=seq, dot_bracket=db)

    def test_find_residue(self):
        ss = secondary_structure.Structure(sequence="AGCU+AGCU",
                                           dot_bracket="((((+))))")
        r = ss.get_residue(1, "A")
        if r is None:
            self.fail("did not find a known residue, find_residue not working")

        r2 = ss.get_residue(uuid=r.uuid)
        if r2 is None:
            self.fail("did not find a known residue, find_residue not working")

        with self.assertRaises(exceptions.SecondaryStructureException):
            ss.get_residue(chain_id="A")

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