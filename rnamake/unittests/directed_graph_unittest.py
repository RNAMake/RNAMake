import unittest
import numpy as np

from rnamake import motif_directed_graph, directed_graph
from rnamake import resource_manager
from rnamake import exceptions

def get_minittr_mdg(rm):
    mdg = motif_directed_graph.MotifDirectedGraph(rm)
    mdg.add_motif_by_name("HELIX.IDEAL.3")
    mdg.add_motif_by_name("GAAA_tetraloop", "A229-A245", 0, 1)
    mdg.add_motif_by_name("HELIX.IDEAL.3", None, 1, 1)
    mdg.add_motif_by_name("HELIX.IDEAL.3", None, 1, 2)
    mdg.add_motif_by_name("HELIX.IDEAL.4", None, 2, 1)
    mdg.add_motif_by_name("TWOWAY.2HW8.0", "B10-B29", 4, 1)
    mdg.add_motif_by_name("HELIX.IDEAL.5", None, 5, 1)
    mdg.add_motif_by_name("TWOWAY.1S72.34", "01529-01662", 6, 1)
    mdg.add_motif_by_name("HELIX.IDEAL.6", None, 7, 1)
    mdg.add_motif_by_name("TWOWAY.1S72.29", "01312-01342", 8, 1)
    mdg.add_motif_by_name("HELIX.IDEAL.5", None, 9, 1)
    mdg.add_motif_by_name("TWOWAY.1S72.100", "930-950", 10, 1)
    return mdg


class DirectedGraphUnittest(unittest.TestCase):

    def test_creation(self):
        dg = directed_graph.DirectedGraph()
        dg.add_node(99, 2)
        dg.add_node(97, 2, parent_index=0)
        dg.add_node(98, 2, parent_index=0, parent_edge_index=0)
        dg.add_node(96, 2)

        self.failUnless(dg.get_parent_index(0) is None)
        self.failUnless(dg.get_parent_index(1) == 0)
        self.failUnless(dg.get_parent_index(2) == 0)
        self.failUnless(dg.get_parent(1) == 99)

        self.failUnless(dg.get_roots() == [0, 3])

        self.failUnless(dg.are_nodes_connected(0, 1) == True)
        self.failUnless(dg.are_nodes_connected(0, 4) == False)

        dg.remove_node(0)
        self.failUnless(dg.get_parent_index(1) is None)

        dg.add_node(95, 2, parent_index=1)
        dg.remove_edge(1, 4, 1, 0)

        # edge no longer exists cannot remove
        with self.assertRaises(ValueError):
            dg.remove_edge(1, 4, 1, 0)

    def test_iter(self):
        dg = directed_graph.DirectedGraph()
        dg.add_node(99, 3)
        dg.add_node(98, 2, parent_index=0, parent_edge_index=1)
        dg.add_node(97, 2)
        dg.add_node(96, 2, parent_index=1, parent_edge_index=1)
        dg.add_node(95, 2, parent_index=2, parent_edge_index=1)
        dg.add_edge(4, 0, 1, 0)

        nis = []
        for ni in dg:
            nis.append(ni)

        self.failUnless(nis == [0, 1, 3, 2, 4])

    def test_add_graph(self):
        dg = directed_graph.DirectedGraph()
        dg.add_node(99, 2)

        dg1 = directed_graph.DirectedGraph()
        dg1.add_node(100, 99)
        dg1.add_node(97, 99, parent_index=0)
        dg1.add_node(98, 99, parent_index=1)
        dg1.add_node(96, 99)
        dg1.add_edge(1, 3, 2, 2)
        dg1.add_edge(0, 3, 3, 3)

        dg.add_graph(dg1, 0, 0, 0, 0)
        self.failUnless(dg.get_parent_index(1) == 0)
        self.failUnless(dg.get_parent_index(4) is None)
        self.failUnless(dg.are_nodes_connected(1, 4))


class MotifDirectedGraphUnittest(unittest.TestCase):

    def setUp(self):
        self._rm = resource_manager.ResourceManager()
        self._m1 = self._rm.get_motif(name="HELIX.IDEAL.2")
        self._m2 = self._rm.get_motif(name="HELIX.IDEAL.2")
        self._m3 = self._rm.get_motif(name="HELIX.IDEAL.2")

    def test_add_motif(self):
        mdg = motif_directed_graph.MotifDirectedGraph(self._rm)
        mdg.add_motif(self._m1)
        mdg.add_motif(self._m2, 0, 1)

        # cannot add same motif again
        with self.assertRaises(ValueError):
            mdg.add_motif(self._m1, 1, 1)

        # parent position does not exist
        with self.assertRaises(exceptions.GraphIndexException):
            mdg.add_motif(self._m3, 2, 2)

        # cannot add motif to end that doesnt exist
        with self.assertRaises(ValueError):
            mdg.add_motif(self._m3, 1, 2)

        # cannot add motif block end add positon
        with self.assertRaises(ValueError):
            mdg.add_motif(self._m3, 1, 0)

        # cannot add motif block end add positon
        with self.assertRaises(ValueError):
            mdg.add_motif(self._m3, 1, parent_end_name="A1-A8")

        # cannot add motif to end that doesnt exist
        with self.assertRaises(ValueError):
            mdg.add_motif(self._m3, 1, parent_end_name="FAKE")

    def test_add_motif_by_name(self):
        mdg = motif_directed_graph.MotifDirectedGraph(self._rm)
        mdg.add_motif_by_name("HELIX.IDEAL.2")
        mdg.add_motif_by_name("HELIX.IDEAL.2", parent_index=0,
                              parent_end_index=1)

        # incorrect motif name
        with self.assertRaises(ValueError):
            mdg.add_motif_by_name("FAKE.MOTIF", parent_index=1,
                                  parent_end_index=1)

        # incorrect motif name
        with self.assertRaises(ValueError):
            mdg.add_motif_by_name("HELIX.IDEAL.2", m_end_name="FAKE",
                                  parent_index=1, parent_end_index=1)

        mdg.add_motif_by_name("HELIX.IDEAL.2", "A1-A8", 1, 1)

        self.failUnless(len(mdg) == 3)

    def test_remove_motif(self):
        mdg = motif_directed_graph.MotifDirectedGraph(self._rm)
        mdg.add_motif(self._m1)
        mdg.add_motif(self._m2, 0, 1)
        mdg.add_motif(self._m3, 1, 1)

        mdg.remove_motif(1)
        self.failUnless(len(mdg) == 2)

        mdg.add_motif(self._rm.get_motif(name="HELIX.IDEAL.2"), 0, 1)
        self.failUnless(len(mdg) == 3)

    def test_merge(self):
        mdg = motif_directed_graph.MotifDirectedGraph(self._rm)
        mdg.add_motif_by_name("HELIX.IDEAL")
        mdg.add_motif_by_name("HELIX.IDEAL", None, 0, 1)
        rs = mdg.get_merged_rna_structure()

        self.failUnless(rs.num_chains() == 2)
        self.failUnless(rs.num_res() == 6)

    def test_add_connection(self):
        mdg = motif_directed_graph.MotifDirectedGraph(self._rm)
        nway = self._rm.get_motif(name="NWAY.1GID.0")
        mdg.add_motif(self._m1)
        mdg.add_motif(nway, 0, 1)
        mdg.add_motif(self._m2, 1, 1)

        ei1 = mdg.get_motif(1).get_end_index("A138-A180")
        ei2 = mdg.get_motif(1).get_end_index("A141-A162")

        # try connecting through 0th end position
        with self.assertRaises(ValueError):
            mdg.add_connection(1, 2, ei1, 1)

        # try connecting thru an already used end position
        with self.assertRaises(ValueError):
            mdg.add_connection(1, 2, ei2, 1)

        self.failUnless(mdg.is_valid_connection(1, 2, 2, 1) == 0)

    def test_replace_sequence(self):
        mdg = motif_directed_graph.MotifDirectedGraph(self._rm)
        mdg.add_motif_by_name("HELIX.IDEAL.15")
        mdg.add_motif_by_name("TWOWAY.3R1C.32", "g4-h5", 0, 1)
        mdg.add_motif_by_name("HELIX.IDEAL.14", None, 1, 1)
        mdg2 = mdg.get_graph_wo_ideal_helices()

        seq = "ACUGAGGAACGUACGACGGCUUACCGUUUAACUA&UAGUUAAACGGUAAGCGGUCGUACGUUCCUCAGU"
        mdg2.update_sequence(seq)

        #dss = mg.secondary_structure()
        #rna_struc = mg.get_structure()
        #mg.to_pdb("org.pdb", renumber=1, close_chain=1)

        # print rna_struc.get_chain(0).get_residue(0).name
        # print dss.get_chain(0).get_residue(0).name

    def test_replace_ideal_helices(self):
        mdg = motif_directed_graph.MotifDirectedGraph(self._rm)
        mdg.add_motif(self._m1)
        mdg.add_motif(self._m2, 0, 1)
        mdg.add_motif(self._m3, 1, 1)

        mdg2 = mdg.get_graph_wo_ideal_helices()

        diff = mdg.last_motif.get_end(1).diff(mdg2.last_motif.get_end(1))
        self.failUnless(diff < 0.75)

    def test_replace_ideal_helices_2(self):
        mdg = motif_directed_graph.MotifDirectedGraph(self._rm)
        m = self._rm.get_motif(name='HELIX.IDEAL.6')
        m.move(np.array([40, 0, 0]))
        mdg.add_motif(m)

        mdg2 = mdg.get_graph_wo_ideal_helices()
        diff = mdg.last_motif.get_end(1).diff(mdg2.last_motif.get_end(1))
        self.failUnless(diff < 0.75)

    def test_copy(self):
        mdg = motif_directed_graph.MotifDirectedGraph(self._rm)
        mdg.add_motif(self._m1)
        mdg.add_motif(self._m2, 0, 1)
        mdg.add_motif(self._m3, 1, 1)

        mdg_copy = motif_directed_graph.MotifDirectedGraph.copy(mdg)
        self.failUnless(len(mdg_copy) == 3)
        self.failUnless(mdg_copy.are_motif_connected(0, 1))

        diff = mdg.last_motif.get_end(1).diff(mdg_copy.last_motif.get_end(1))
        self.failUnless(diff < 0.001)

    def test_minittr(self):
        mdg = get_minittr_mdg(self._rm)
        diff = mdg.get_motif(1).get_end(2).diff(mdg.last_motif.get_end(1))
        self.failUnless(diff < 27)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
