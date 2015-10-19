import unittest
import rnamake.motif_graph as motif_graph
import rnamake.graph as graph
import rnamake.util as util
import rnamake.eternabot.sequence_designer as sd
import build


class MotifGraphUnittest(unittest.TestCase):

    def test_creation(self):
        mg = motif_graph.MotifGraph()

    def test_add(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)

        mg = motif_graph.MotifGraph()
        for n in mt.graph.nodes:
            mg.add_motif(n.data)

        if len(mg.graph) != 3:
            self.fail("did not get the right number of motifs")

    def test_remove(self):
        builder = build.BuildMotifTree()
        mt = builder.build(5)
        mg = motif_graph.MotifGraph()

        for n in mt.graph.nodes:
            mg.add_motif(n.data)

        mg.write_pdbs()
        mg.remove_motif(0)
        #print len(mg.structure.ends)
        #mg.structure.to_pdb("test.pdb")

    def test_copy(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.graph.nodes:
            mg.add_motif(n.data)

        new_mg = mg.copy()

        if len(mg.graph) != len(new_mg.graph):
            self.fail("did not copy correctly")

        new_mg.remove_motif(1)

    def test_replace_ideal_helices(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.graph.nodes:
            mg.add_motif(n.data)

        new_mg = mg.copy()
        new_mg.replace_ideal_helices()

        d1 = mg.last_node().data.ends[1].d()
        d2 = new_mg.last_node().data.ends[1].d()
        diff = util.distance(d1, d2)
        if diff > 1:
            self.fail("replacing ideal helices messed up graph")

    def _test_get_ss(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.graph.nodes:
            mg.add_motif(n.data)

        new_mg = mg.copy()
        new_mg.replace_ideal_helices()
        ss = new_mg.designable_secondary_structure(leaf_pos=3)
        designer = sd.SequenceDesigner()
        r = designer.design(ss.dot_bracket(), ss.sequence())
        ss.replace_sequence(r[0].sequence)
        new_mg.replace_helix_sequence(ss)
        #new_mg.write_pdbs()

    def test_get_pose(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.graph.nodes:
            mg.add_motif(n.data)

        new_mg = mg.copy()
        new_mg.replace_ideal_helices()
        ss = new_mg.designable_secondary_structure(leaf_pos=3)

    def test_chain_combine(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.graph.nodes:
            mg.add_motif(n.data)

        #mg.structure.to_pdb("test.pdb")
        #ss = mg.designable_secondary_structure()

        #new_mg = mg.copy()
        #new_mg.replace_ideal_helices()

    def test_secondary_structure(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.graph.nodes:
            mg.add_motif(n.data)

        print mg.secondary_structure()

    def test_designable_secondary_structure(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.graph.nodes:
            mg.add_motif(n.data)


        mg.replace_ideal_helices()
        #mg.write_pdbs()
        #for i, c in enumerate(mg.structure.chains()):
        #    c.to_pdb("c."+str(i)+".pdb")

        #mg.structure.to_pdb("test.pdb")

        ss = mg.designable_secondary_structure()

        designer = sd.SequenceDesigner()
        r = designer.design(ss.dot_bracket(), ss.sequence())
        ss.replace_sequence(r[0].sequence)
        print ss
        mg.replace_helix_sequence(ss)
        mg.write_pdbs()
        mg.structure.to_pdb("test.pdb")


def main():
    unittest.main()

if __name__ == '__main__':
    main()