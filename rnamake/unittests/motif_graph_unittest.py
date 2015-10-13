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

    def test_copy(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.graph.nodes:
            mg.add_motif(n.data)

        new_mg = mg.copy()
        if len(mg.graph) != len(new_mg.graph):
            self.fail("did not copy correctly")

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

    def test_get_ss(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.graph.nodes:
            mg.add_motif(n.data)

        new_mg = mg.copy()
        new_mg.replace_ideal_helices()
        ss = new_mg.designable_secondary_structure()
        designer = sd.SequenceDesigner()
        r = designer.design(ss.dot_bracket(), ss.sequence())
        ss.replace_sequence(r[0].sequence)
        new_mg.replace_helix_sequence(ss)




def main():
    unittest.main()

if __name__ == '__main__':
    main()