import tree
import secondary_structure_factory as ssf
import secondary_structure


class MotifTreeTopology(object):
    def __init__(self, connectivity):
        self.tree = tree.TreeDynamic()

        offset = 0
        for i, c in enumerate(connectivity):
            if c[3] != "":
                if c[3][0:5] != "HELIX":
                    parent_end_ss_id = c[1]
                    if i > 0:
                        n = self.tree.get_node(int(c[2])+offset)
                        if n.data.motif_name[2] == "=":
                            ss = ssf.ss_id_to_secondary_structure(n.data.end_ss_id)
                            parent_end_ss_id = secondary_structure.assign_end_id(ss, ss.ends[1])
                    self.tree.add_data(MotifTreeTopologyNodeData(c[3], c[0], parent_end_ss_id),
                                       int(c[2])+offset)
                else:
                    ss = ssf.ss_id_to_secondary_structure(c[0])
                    conn = ss.motif_topology_from_end(ss.ends[0])
                    start = offset+int(c[2])
                    for j, c2 in enumerate(conn):
                        spl = c2[0].split("_")
                        name = spl[0][0]+spl[2][1]+"="+spl[0][1]+spl[2][0]
                        if j == 0:
                            d = MotifTreeTopologyNodeData(name, c2[0], c[1])

                        else:
                            d = MotifTreeTopologyNodeData(name, c2[0], c2[1])
                        self.tree.add_data(d, start+j)

                    offset += len(conn)-1



    def __repr__(self):
        s = "<(MotifTreeTopology: #nodes: %d\n" % (len(self.tree))
        for n in self.tree.nodes:

            s += "\t index: %d name: %s parent %s\n" % (n.index, n.data.motif_name,
                                                        n.parent_index())
        s += ")"

        return s


class MotifTreeTopologyNodeData(object):
    def __init__(self, motif_name, end_ss_id, parent_end_ss_id):
        self.motif_name, self.end_ss_id = motif_name, end_ss_id
        self.parent_end_ss_id = parent_end_ss_id




