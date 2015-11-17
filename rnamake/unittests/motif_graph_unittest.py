import unittest
import rnamake.motif_graph as motif_graph
import rnamake.motif_type as motif_type
import rnamake.graph as graph
import rnamake.util as util
import rnamake.eternabot.sequence_designer as sd
import rnamake.resource_manager as rm
from rnamake import motif_topology
import build


class MotifGraphUnittest(unittest.TestCase):

    def test_creation(self):
        mg = motif_graph.MotifGraph()

    def test_add(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)

        mg = motif_graph.MotifGraph()
        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        if len(mg.graph) != 3:
            self.fail("did not get the right number of motifs")

        mg.write_pdbs()
        mg.merger.to_pdb("test.pdb")

    def test_remove(self):
        builder = build.BuildMotifTree()
        mt = builder.build(2)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)


        exit()
        mg.remove_motif(1)
        # mg.write_pdbs()
        #mg.merger.to_pdb("test.pdb")
        #print len(mg.structure.ends)
        #mg.structure.to_pdb("test.pdb")

    def test_copy(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        new_mg = mg.copy()

        if len(mg.graph) != len(new_mg.graph):
            self.fail("did not copy correctly")

        new_mg.remove_motif(1)

    def test_replace_ideal_helices(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        mg.write_pdbs("org")
        new_mg = mg.copy()
        new_mg.replace_ideal_helices()
        new_mg.write_pdbs()
        new_mg.merger.to_pdb("test.pdb")

        d1 = mg.last_node().data.ends[1].d()
        d2 = new_mg.last_node().data.ends[1].d()
        diff = util.distance(d1, d2)
        if diff > 1:
            self.fail("replacing ideal helices messed up graph")

    def _test_get_ss(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        new_mg = mg.copy()
        new_mg.replace_ideal_helices()
        ss = new_mg.designable_secondary_structure()
        designer = sd.SequenceDesigner()
        r = designer.design(ss.dot_bracket(), ss.sequence())
        ss.replace_sequence(r[0].sequence)
        new_mg.replace_helix_sequence(ss)
        #new_mg.write_pdbs()

    def test_get_pose(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        new_mg = mg.copy()
        new_mg.replace_ideal_helices()
        ss = new_mg.designable_secondary_structure()

    def test_secondary_structure(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()
        for n in mt.tree.nodes:
            mg.add_motif(n.data)
        #mg.write_pdbs("org")
        mg.replace_ideal_helices()

        mg.write_pdbs()
        ss = mg.designable_secondary_structure()

        designer = sd.SequenceDesigner()
        r = designer.design(ss.dot_bracket(), ss.sequence())
        ss.replace_sequence(r[0].sequence)
        mg.replace_helix_sequence(ss)
        #mg.write_pdbs("new")
        #mg.merger.to_pdb("test.pdb")
        #mg.write_pdbs()

    def test_designable_secondary_structure(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)


        mg.write_pdbs("org")
        mg.replace_ideal_helices()
        #mg.write_pdbs()
        #for i, c in enumerate(mg.structure.chains()):
        #    c.to_pdb("c."+str(i)+".pdb")

        #mg.structure.to_pdb("test.pdb")
        mg.write_pdbs()
        ss = mg.designable_secondary_structure()

        designer = sd.SequenceDesigner()
        r = designer.design(ss.dot_bracket(), ss.sequence())
        ss.replace_sequence(r[0].sequence)
        mg.replace_helix_sequence(ss)
        mg.write_pdbs("new")
        mg.merger.to_pdb("test.pdb")

    def test_designable_secondary_structure_2(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        hairpin = rm.manager.mlibs['hairpin'].get_random()
        mg.add_motif(hairpin)
        mg.write_pdbs()

    def test_ss_real_case(self):
        rm.manager.add_motif("resources/motifs/tetraloop_receptor_min")
        mg = motif_graph.MotifGraph()
        mg.add_motif(rm.manager.get_motif(name="tetraloop_receptor_min",
                                          end_name="A228-A246"))
        mg.add_motif(rm.manager.get_motif(name="HELIX.IDEAL.20"), parent_end_name="A221-A252")
        mg.add_motif(rm.manager.get_motif(name="tetraloop_receptor_min",
                                          end_name="A221-A252"))

        mg.replace_ideal_helices()
        ss = mg.designable_secondary_structure()
        build.fill_basepairs_in_ss(ss)
        mg.replace_helix_sequence(ss)
        #print ss
        #mg.write_pdbs()

    def test_to_tree(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        #mg.write_pdbs("org")
        #mg.secondary_structure()

        mt2 = motif_topology.graph_to_tree(mg)
        mt2.write_pdbs()

def main():
    unittest.main()

if __name__ == '__main__':
    main()














