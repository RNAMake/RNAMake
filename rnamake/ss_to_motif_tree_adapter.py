import motif_tree_topology
import motif_tree
import graph
import resource_manager

class SStoMotifTreeAdapter(object):
    def __init__(self):
        pass

    def convert(self, sstree):
        mtt = motif_tree_topology.MotifTreeTopology(sstree)

        for n in mtt:
            print n.data.ends[0].ss_data

        mt = motif_tree.MotifTree()

        i = -1
        seen_nodes = {}
        for n in graph.transverse_graph(mtt.graph, 0):
            i += 1
            if i == 0:
                avail_pos = n.available_children_pos()
                if len(avail_pos) == 0:
                    raise ValueError("cannot start build mt from mtt from this node")
                m = resource_manager.manager.get_motif(n.data.ends[ avail_pos[0] ].get_id())
                index = mt.add_motif(m)
                seen_nodes[n.index] = index


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
                avail_pos = mt.nodes[mt_index].available_indexes()
                for i in range(len(mt.nodes[mt_index].motif.ends)):
                    if mt.nodes[mt_index].motif.get_end_id(i) == correct_id and i in avail_pos:
                        mt_parent_index = i
                        break

                m = resource_manager.manager.get_motif(n.data.ends[ end_index ].get_id())
                print m.name, mt_parent_index, parent_end_index, 0
                parent = mt.nodes[mt_index]
                mt_n = mt.add_motif(m, parent=parent, parent_index=mt_parent_index, end_index=0, end_flip=0)
                seen_nodes[n.index] = mt_n.index
                mt.write_pdbs()


