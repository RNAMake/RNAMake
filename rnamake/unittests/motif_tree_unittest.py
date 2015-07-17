import unittest
import build
import random
import rnamake.motif_tree as motif_tree
import rnamake.resource_manager as rm
import rnamake.secondary_structure_factory as ssfactory
import rnamake.motif_type as motif_type
import rnamake.eternabot.sequence_designer as sequence_designer
import rnamake.motif_tree_topology as motif_tree_topology
import util

class MotifTreeUnittest(unittest.TestCase):

    def test_creation(self):
        mt = rnamake.motif_tree.MotifTree()

    def test_add_motif(self):
        mt = rnamake.motif_tree.MotifTree()
        m = rm.manager.get_motif("HELIX.IDEAL.2")
        mt.add_motif(m)
        mt.add_motif(m)
        if len(mt) != 2:
            self.fail("did not add motifs properly")

    def test_remove_node(self):
        mt = rnamake.motif_tree.MotifTree()
        m = rm.manager.get_motif("HELIX.IDEAL.2")
        mt.add_motif(m)
        mt.add_motif(m)

        mt.remove_node()
        if len(mt) != 1:
            self.fail("did not remove node correctly")

    def test_merge(self):
        mt = rnamake.motif_tree.MotifTree()
        m = rm.manager.get_motif("HELIX.IDEAL.2")
        mt.add_motif(m)
        mt.add_motif(m)
        p = mt.to_pose(chain_closure=0)


        """return
        f = open("test.out")
        s = f.readline()
        f.close()

        mt = rnamake.motif_tree.str_to_motif_tree(s)
        pose = mt.get_pose()
        if len(pose.ends) != 2:
            self.fail("did not merge properly")"""

    def test_readd(self):
        mt = rnamake.motif_tree.MotifTree()
        rm = rnamake.resource_manager.ResourceManager()
        m = rm.get_motif("HELIX.IDEAL")
        mt.add_motif(m)
        m = mt.nodes[1].motif
        mt.remove_node(mt.last_node)
        node = mt.add_motif(m)
        print node
        mt.write_pdbs()

    def test_replace_ideal_helices(self):
        builder = build.BuildMotifTree()
        mt = builder.build()
        mt.write_pdbs()
        p = mt.to_pose()
        seq = p.optimized_sequence()
        db = p.dot_bracket()

        ss = ssfactory.factory.get_structure(seq, db)
        connectivity = ss.motif_topology_from_end(ss.ends[0])
        mt1 = rnamake.motif_tree.motif_tree_from_topology(connectivity)

    def test_secondary_structure(self):
        rm.manager.add_motif("resources/motifs/tetraloop_receptor_min")
        mt = motif_tree.MotifTree()
        mt.add_motif(rm.manager.get_motif(name="tetraloop_receptor_min",
                                          end_name="A228-A246"))
        mt.add_motif(rm.manager.get_motif(name="HELIX.IDEAL.20"), parent_end_name="A221-A252")
        mt.add_motif(rm.manager.get_motif(name="tetraloop_receptor_min",
                                          end_name="A228-A246"))
        mt.write_pdbs("org")
        p = mt.to_pose()
        ss = p.designable_secondary_structure()
        for ss_r in ss.residues():
            r = p.get_residue(uuid=ss_r.uuid)
            if r is None:
                self.fail("did nto copy secondary structure properly")

        pairs = ["AU", "UA", "GC", "CG"]
        for bp in ss.basepairs:
            if bp.res1.name == "N" and bp.res2.name == "N":
                p = random.choice(pairs)
                bp.res1.name = p[0]
                bp.res2.name = p[1]

        conn =  ss.motif_topology_from_end(ss.ends[0])
        mtt = motif_tree_topology.MotifTreeTopology(conn)
        #for n in mtt.tree:
        #    print n.data.motif_name, n.data.parent_end_ss_id
        mt2 = motif_tree.motif_tree_from_topology_2(mtt)
        mt2.write_pdbs()







def main():
    unittest.main()

if __name__ == '__main__':
    main()
