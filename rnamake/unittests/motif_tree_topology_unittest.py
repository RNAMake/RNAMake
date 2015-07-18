import unittest
import random
import rnamake.sqlite_library as sqlite_library
import rnamake.motif_tree_topology as motif_tree_topology
import rnamake.motif_tree as motif_tree
import rnamake.resource_manager as resource_manager
import rnamake.motif_type
import rnamake.settings as settings
import rnamake.motif_factory as motif_factory
import rnamake.secondary_structure as secondary_structure
import rnamake.ss_tree as ss_tree
import build
import rnamake.motif_type as motif_type

class MotifTreeTopologyUnittest(unittest.TestCase):

    def _fill_basepairs_in_ss(self, ss):
        pairs = ["AU", "UA", "GC", "CG"]
        for bp in ss.basepairs:
            if bp.res1.name == "N" and bp.res2.name == "N":
                p = random.choice(pairs)
                bp.res1.name = p[0]
                bp.res2.name = p[1]

    def test_creation(self):
        builder = build.BuildSecondaryStructure()
        ss = builder.build_helix()
        conn = ss.motif_topology_from_end(ss.ends[1])
        mtt = motif_tree_topology.MotifTreeTopology(conn)
        mt = motif_tree.motif_tree_from_topology_2(mtt)
        if len(mt) != len(ss.motifs('ALL')):
            self.fail("did not build properly")

    def test_twoway(self):
        builder = build.BuildMotifTree()
        mt = builder.build(size=10)
        ss = mt.designable_secondary_structure()
        self._fill_basepairs_in_ss(ss)
        con = ss.motif_topology_from_end(ss.ends[1])
        mtt = motif_tree_topology.MotifTreeTopology(con)
        mt2 = motif_tree.motif_tree_from_topology_2(mtt)

    def test_nway(self):
        builder = build.BuildMotifTree(lib_names=['ideal_helices', 'nway'])
        mt = builder.build_specific('HELIX.IDEAL.3,NWAY.2VQE.9,HELIX.IDEAL.3'.split(','))
        #mt = builder.build(size=3)
        for n in mt:
            print n.data.name
        mt.write_pdbs("org")
        ss = mt.designable_secondary_structure()
        self._fill_basepairs_in_ss(ss)
        con = ss.motif_topology_from_end(ss.ends[1])
        mtt = motif_tree_topology.MotifTreeTopology(con)
        #for n in mtt.tree.nodes:
        #    print n.index, n.data.motif_name, n.data.parent_end_ss_id, n.parent_index()
        mt2 = motif_tree.motif_tree_from_topology_2(mtt)
        mt2.write_pdbs()

    def _test_build_mt(self):
        sstree = ss_tree.SS_Tree("(((+)))", "GAG+UUC")
        mtt = motif_tree_topology.MotifTreeTopology(sstree)
        mt = motif_tree.MotifTree()

        i = -1
        seen_nodes = {}
        for n in graph.transverse_graph(mtt.graph, 1):
            i += 1
            if i == 0:
                avail_pos = n.available_children_pos()
                if len(avail_pos) == 0:
                    raise ValueError("cannot start build mt from mtt from this node")
                m = resource_manager.manager.get_motif(n.data.ends[ avail_pos[0] ].get_id())
                mt_n = mt.add_motif(m, end_index=0, end_flip=0)
                seen_nodes[n.index] = mt_n.index

            else:
                parent_index = -1
                mt_index = -1
                for pi, mt_i in seen_nodes.iteritems():
                    if n.is_connected(pi):
                        if mt_i > mt_index:
                            mt_index = mt_i
                            parent_index = pi

                end_index, parent_end_index = n.end_indexes_with_node(parent_index)
                parent = mtt.graph.get_node(parent_index)
                correct_id = parent.data.ends[ parent_end_index ].get_id()

                mt_parent_index = -1
                for i in range(len(mt.nodes[mt_index].motif.ends)):
                    if mt.nodes[mt_index].motif.get_end_id(i) == correct_id:
                        mt_parent_index = i
                        break


                m = resource_manager.manager.get_motif(n.data.ends[ end_index ].get_id())
                parent = mt.nodes[mt_index]
                mt.add_motif(m, parent=parent, end_index=0, parent_index=mt_parent_index, end_flip=0)

        print len(mt.nodes)
        mt.write_pdbs()


def main():
    unittest.main()

if __name__ == '__main__':
    main()




















