import unittest
import build

from rnamake import motif_topology, settings, motif_graph
from rnamake import resource_manager as rm


class MotifTopologyUnittest(unittest.TestCase):

    def test_graph_to_tree(self):
        builder = build.BuildMotifGraph()
        mg = builder.build(3)

        c = motif_topology.GraphtoTree()
        mt = c.convert(mg)
        self.failIf(len(mt) != len(mg))

        c = motif_topology.GraphtoTree()
        mt = c.convert(mg, start=mg.last_node(), start_end_index=1)

        self.failIf(len(mt) != len(mg))

        mg.replace_ideal_helices()
        dss = mg.designable_secondary_structure()
        build.fill_basepairs_in_ss(dss)
        mg.replace_helix_sequence(dss)
        #print dss.sequence()
        #print mg.sequence()

        c = motif_topology.GraphtoTree()
        bpoints = mg.get_build_points()

        mt = c.convert(mg, start=bpoints[0].node,
                       start_end_index=bpoints[0].end_index)

        self.failIf(len(mt) != len(mg))

    def test_last_end(self):
        path = settings.UNITTEST_PATH + "/resources/motif_graph/mini_ttr.mg"
        f = open(path)
        lines = f.readlines()
        f.close()

        mg = motif_graph.MotifGraph(mg_str=lines[0])
        """
        for n in mg.graph:
            n.data.new_res_uuids()

        c = motif_topology.GraphtoTree()
        mt = c.convert(mg)"""

        n = mg.get_node(m_name = "GAAA_tetraloop")
        last_node = n.connections[1].partner(n.index)

        c = motif_topology.GraphtoTree()
        mt = c.convert(mg, last_node=last_node)
        #mt.write_pdbs()












def main():
    unittest.main()

if __name__ == '__main__':
    main()
