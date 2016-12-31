import unittest
import numpy as np
from rnamake import motif_state_graph, motif_graph, exceptions, util, settings
from rnamake import resource_manager as rm


class MotifStateGraphUnittest(unittest.TestCase):

    def test_create(self):
        msg = motif_state_graph.MotifStateGraph()

    def test_add_state(self):
        msg = motif_state_graph.MotifStateGraph()

        ms1 = rm.manager.get_state(name="HELIX.IDEAL.2")
        ms2 = rm.manager.get_state(name="HELIX.IDEAL.2")
        msg.add_state(ms1)

        # can only add the motif state with the same unique indenitifer twice
        with self.assertRaises(exceptions.MotifStateGraphException):
            msg.add_state(ms1)

        # can never use parent_end_index=0 for a graph as that is where that node
        # is already connected to another node
        with self.assertRaises(exceptions.MotifStateGraphException):
            msg.add_state(ms2, parent_end_index=0)

        # must supply a motif or motif name
        with self.assertRaises(exceptions.MotifStateGraphException):
            msg.add_state()

        # motif not found in resource manager
        with self.assertRaises(exceptions.MotifStateGraphException):
            msg.add_state(m_name="FAKE")

        # catches invalid parent_index
        with self.assertRaises(exceptions.MotifStateGraphException):
            msg.add_state(ms2, parent_index=2)

        # invalid parent_end_index, has only 0 and 1
        with self.assertRaises(exceptions.MotifStateGraphException):
            msg.add_state(ms2, parent_end_index=3)

        # invalid parent_end_name, is the name of end 0
        with self.assertRaises(exceptions.MotifStateGraphException):
            msg.add_state(ms2, parent_end_name="A4-A5")

        # invalid parent_end_name, cannot be found as an end in motif
        with self.assertRaises(exceptions.MotifStateGraphException):
            msg.add_state(ms2, parent_end_name="FAKE")

    def test_setup_from_mg(self):
        mg = motif_graph.MotifGraph()
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m3 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        nway = rm.manager.get_motif(name="NWAY.1GID.0")
        mg.add_motif(m1)
        mg.add_motif(nway)
        mg.add_motif(m2)
        mg.add_motif(m3, parent_index=1)
        mg.add_connection(2, 3)

        msg = motif_state_graph.MotifStateGraph(mg)
        self.failUnless(len(msg) == len(mg))
        self.failUnless(msg.get_node(1).connections[2].partner(1).index == 3)
        self.failUnless(len(msg.graph.connections) == 4)

    def test_setup_from_mg_2(self):
        mg = motif_graph.MotifGraph()
        for i in range(10):
            mg.add_motif(m_name="HELIX.IDEAL")

        msg = motif_state_graph.MotifStateGraph(mg)

        d1 = mg.get_node(9).data.ends[1].d()
        d2 = msg.get_node(9).data.cur_state.end_states[1].d
        self.failUnless(util.distance(d1, d2) < 0.1)

    def test_to_motif_graph(self):
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        tc = rm.manager.get_motif(name="TC.1S72.0")

        mg = motif_graph.MotifGraph()
        mg.add_motif(tc)
        mg.add_motif(m1)
        mg.add_motif(m2, parent_index=0)
        msg = motif_state_graph.MotifStateGraph(mg)
        mg2 = msg.to_motif_graph()

        self.failUnless(len(mg2) == len(mg))
        atoms1 = mg.get_structure().structure.atoms()
        atoms2 = mg2.get_structure().structure.atoms()

        for i in range(len(atoms1)):
            diff = util.distance(atoms1[i].coords, atoms2[i].coords)
            if diff > 0.1:
                print atoms1[i], atoms2[i]
                self.fail("atoms are not equal")


        mg.add_connection(0, 2)
        msg = motif_state_graph.MotifStateGraph(mg)
        mg2 = msg.to_motif_graph()

        self.failUnless(len(mg.get_structure().residues()) == \
                        len(mg2.get_structure().residues()))

    def test_to_motif_graph_2(self):
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2.move(np.array([40, 0, 0]))

        mg = motif_graph.MotifGraph()
        mg.add_motif(m1)
        mg.add_motif(m2, orphan=1)
        mg.add_motif(m_name="HELIX.IDEAL.2")
        mg.add_motif(m_name="HELIX.IDEAL.2")
        mg.add_motif(m_name="HELIX.IDEAL.2")
        mg.add_motif(m_name="HELIX.IDEAL.2", parent_index=0)
        mg.add_motif(m_name="HELIX.IDEAL.2")
        mg.add_motif(m_name="HELIX.IDEAL.2")

        msg = motif_state_graph.MotifStateGraph(mg)
        mg2 = msg.to_motif_graph()

        self.failUnless(len(mg) == len(mg2))

        mg.add_connection(4, 7)
        msg = motif_state_graph.MotifStateGraph(mg)
        mg2 = msg.to_motif_graph()

        self.failUnless(len(mg2.get_structure().chains()) == 2)

    def test_add_connection(self):
        m1 = rm.manager.get_state(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_state(name="HELIX.IDEAL.2")
        tc = rm.manager.get_motif(name="TC.1S72.0").get_state()

        msg = motif_state_graph.MotifStateGraph()
        msg.add_state(tc)
        msg.add_state(m1)
        msg.add_state(m2, parent_index=0)
        msg.add_connection(0, 2)

        mg = msg.to_motif_graph()
        self.failUnless(len(mg.get_structure().chains()) == 1)

    def test_remove_state(self):
        ms1 = rm.manager.get_state(name="HELIX.IDEAL.2")
        ms2 = rm.manager.get_state(name="HELIX.IDEAL.2")
        ms3 = rm.manager.get_state(name="HELIX.IDEAL.2")
        msg = motif_state_graph.MotifStateGraph()

        msg.add_state(ms1)
        msg.add_state(ms2)
        msg.add_state(ms3)

        msg.remove_state(1)

        self.failUnless(len(msg) == 2)
        msg.remove_state(2)
        msg.add_state(ms2)
        msg.add_state(ms3)

        mg = msg.to_motif_graph()
        self.failUnless(len(mg.get_structure().chains()) == 2)

    def test_remove_node_level(self):
        ms1 = rm.manager.get_state(name="HELIX.IDEAL.2")
        ms2 = rm.manager.get_state(name="HELIX.IDEAL.2")
        ms3 = rm.manager.get_state(name="HELIX.IDEAL.2")
        msg = motif_state_graph.MotifStateGraph()

        msg.add_state(ms1)
        msg.increase_level()
        msg.add_state(ms2)
        msg.add_state(ms3)
        msg.remove_node_level(1)
        self.failUnless(len(msg) == 1)

    def test_replace_state(self):
        ms1 = rm.manager.get_state(name="HELIX.IDEAL.2")
        ms2 = rm.manager.get_state(name="TWOWAY.2PN4.4")
        ms3 = rm.manager.get_state(name="HELIX.IDEAL.2")
        msg = motif_state_graph.MotifStateGraph()

        msg.add_state(ms1)
        msg.add_state(ms2)
        msg.add_state(ms3)

        msg.replace_state(1, rm.manager.get_state(name="HELIX.IDEAL.2"))
        self.failUnless(msg.get_node(1).data.name() == "HELIX.IDEAL.2")

    def test_multiple_alignments(self):
        base_dir = settings.UNITTEST_PATH + "resources/motif_graph/"
        f = open(base_dir+"tecto_chip_only.mg")
        l = f.readline()
        f.close()

        mg = motif_graph.MotifGraph(mg_str=l)
        msg = motif_state_graph.MotifStateGraph(mg)
        mg2 = msg.to_motif_graph()

        self.failUnless(mg.sequence() == mg2.sequence())


def main():
    unittest.main()

if __name__ == '__main__':
    main()