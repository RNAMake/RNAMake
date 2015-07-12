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

class SecondaryStructureUnittest(unittest.TestCase):

    def _test_assign_secondary_structure(self):
        builder = build.BuildMotifTree()
        mt = builder.build()
        #for n in mt:
        #    print n.data.name

        p = mt.to_pose()
        #print p.secondary_structure()
        #print p.sequence()
        #mt.write_pdbs()

    def test_creation(self):
        ss = secondary_structure.SecondaryStructure("AGCU+AGCU","((((+))))")

    def test_find_residue(self):
        ss = secondary_structure.SecondaryStructure("AGCU+AGCU","((((+))))")
        r = ss.get_residue(0, "A")
        if r is None:
            self.fail("did not find a known residue, finde_residue not working")

    def test_to_str(self):
        builder = build.BuildSecondaryStructure()
        ss = builder.build_helix(10)
        s = ss.to_str()
        ss1 = secondary_structure.str_to_secondary_structure(s)

        if len(ss.residues()) != len(ss1.residues()):
            self.fail("did not recover all residues from str_to_secondary_structure")

        if len(ss.elements["BP_STEP"]) != len(ss1.elements["BP_STEP"]):
            self.fail("did not get all basepair steps again")

    def test_copy(self):
        builder = build.BuildSecondaryStructure()
        ss = builder.build_helix(2)

        ss_copy = ss.copy()
        ss_copy = ss_copy.copy()

    def test_parse(self):
        seq = "UG&CA&CGACACAG"
        db  = "((&))&(......)"
        ss = ssfactory.factory.get_structure(seq, db)
        for bp in ss.basepairs:
            print bp.res1.num, bp.res2.num

    def test_motif_topology_from_end(self):
        builder = build.BuildSecondaryStructure()
        ss = builder.build_helix(10)
        connectivity = ss.motif_topology_from_end(ss.ends[0])
        for c in connectivity:
            print c
        mt = motif_tree.motif_tree_from_topology(connectivity)
        mt.write_pdbs()

def main():
    unittest.main()

if __name__ == '__main__':
    main()