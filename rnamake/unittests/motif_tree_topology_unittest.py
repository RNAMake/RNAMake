import unittest
import random
import rnamake.graph as graph
import rnamake.ss_tree as ss_tree
import rnamake.motif_tree_topology as motif_tree_topology
import rnamake.motif_tree as motif_tree
import rnamake.resource_manager as resource_manager
import rnamake.motif_library_sqlite
import rnamake.motif_type
import rnamake.settings as settings
import rnamake.ss_to_motif_tree_adapter as ss_to_motif_tree_adapter

class MotifTreeTopologyUnittest(unittest.TestCase):

    def _build_random_mt(self):
        h_lib = rnamake.motif_library_sqlite.MotifLibrarySqlite(rnamake.motif_type.HELIX)
        path = settings.RESOURCES_PATH + "motif_libraries/twoway_aligned.db"
        twoways = rnamake.motif_library_sqlite.MotifLibrarySqlite(libpath=path)
        twoways.mtype = rnamake.motif_type.TWOWAY
        twoways.load_all()

        mt = motif_tree.MotifTree()
        for i in range(2):
            if i % 2 == 0:
                num = random.randint(1,20)
                m = h_lib.get_motif("HELIX.IDEAL."+str(num))
                mt.add_motif(m, end_index=1, end_flip=0)
            else:
                m = twoways.random_motif()
                mt.add_motif(m, end_index=m.end_to_add, end_flip=0)

        return mt

    def _specific_mt_build(self):
        mt = motif_tree.MotifTree()
        m = resource_manager.manager.get_motif("HELIX.IDEAL.2")
        mt.add_motif(m, end_index=1, end_flip=0)
        path = settings.RESOURCES_PATH + "motif_libraries/twoway_aligned.db"
        twoways = rnamake.motif_library_sqlite.MotifLibrarySqlite(libpath=path)
        twoways.mtype = rnamake.motif_type.TWOWAY
        m = twoways.get_motif("TWOWAY.2ZY6.0-0")
        mt.add_motif(m, end_index=m.end_to_add, end_flip=0)
        return mt


    def test_creation(self):
        #sstree = ss_tree.SS_Tree("(((+)))", "GAG+UUC")
        sstree = ss_tree.SS_Tree("((.((+)).((+)).))", "GGAGG+CCAGG+CCACC")

        mtt = motif_tree_topology.MotifTreeTopology(sstree)
        if len(mtt) != 4:
            raise ValueError("did not get the correct number of nodes")

        mt = self._specific_mt_build()
        for n in mt.nodes:
            print n.motif.name

        p = mt.to_pose()

        print p.secondary_structure()
        print p.sequence()

        sstree = ss_tree.SS_Tree(p.secondary_structure(), p.sequence())

        adapter = ss_to_motif_tree_adapter.SStoMotifTreeAdapter()
        mt = adapter.convert(sstree)




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
