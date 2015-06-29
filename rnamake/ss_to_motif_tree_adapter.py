import motif_tree_topology
import motif_tree
import graph
import resource_manager as rm

class SStoMotifTreeAdapter(object):
    def __init__(self):
        pass

    def convert(self, sstree):
        mtt = motif_tree_topology.MotifTreeTopology(sstree)

        #for n in mtt:
        #    print n.data.ends[0].get_id()

        mt = motif_tree.MotifTree()

        i = -1
        seen_nodes = {}
        for n in graph.transverse_graph(mtt.graph, 0):
            i += 1
            if i == 0:
                avail_pos = n.available_children_pos()
                if len(avail_pos) == 0:
                    raise ValueError("cannot start build mt from mtt from this node")
                end_id = n.data.ends [ avail_pos[0] ].get_id()
                m = rm.manager.get_motif(ss_id=end_id)
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
                avail_pos = mt.get_node(mt_index).available_children_pos()
                for p in avail_pos:
                    if mt.get_node(mt_index).data.end_ids[p] == correct_id:
                        mt_parent_index = p
                        break

                m = rm.manager.get_motif(ss_id=n.data.ends[ end_index ].get_id())
                mt_n = mt.add_motif(m, mt_index, parent_end_index)

                seen_nodes[n.index] = mt_n

        return mt


