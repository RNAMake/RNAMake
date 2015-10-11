import tree
import secondary_structure_factory as ssf
import secondary_structure

class ConnectivityNode(object):
    def __init__(self, end_id, parent_end_id, parent_index, name, start_pos):
        self.end_id, self.parent_end_id = end_id, parent_end_id
        self.parent_index, self.name, self.start_pos = parent_index, name, start_pos
        self.size = 0
        self.end_pos = 0

    def copy(self):
        return ConnectivityNode(self.end_id, self.parent_end_id, self.parent_index,
                                self.name, self.start_pos)

    def __repr__(self):
        s = "<(ConnectivityNode: name: %s, s_pos: %s, e_pos: %s, parent_index: %s\n" % \
            (self.name, self.start_pos, self.end_pos, self.parent_index)
        s+= "                    end_id: %s, parent_end_id: %s, size: %s )>" % \
            (self.end_id, self.parent_end_id, self.size)

        return s

class MotifTreeTopology(object):
    def __init__(self, connectivity):
        self.tree = tree.TreeDynamic()

        conn_nodes = [ConnectivityNode(*c) for c in connectivity]
        new_conn_nodes = []

        for i, c in enumerate(conn_nodes):

            new_c = c.copy()

            if new_c.name[0:5] == "HELIX":
                spl = new_c.name.split(".")
                if len(spl) == 3:
                    num = int(spl[2])
                else:
                    num = 0
                new_c.size = num

            if i == 0:
                new_c.end_pos = new_c.start_pos + new_c.size
                new_conn_nodes.append(new_c)
                continue
            parent = new_conn_nodes[c.parent_index]

            new_c.parent_index = parent.end_pos
            new_c.start_pos = new_conn_nodes[-1].start_pos + new_conn_nodes[-1].size + 1
            new_c.end_pos = new_c.start_pos + new_c.size

            new_conn_nodes.append(new_c)

        #for c in new_conn_nodes:
        #    print c
        #exit()

        for i, c in enumerate(new_conn_nodes):
            if c.name != "":
                if c.name[0:5] != "HELIX":
                    parent_end_ss_id = c.parent_end_id
                    if i > 0:
                        n = self.tree.get_node(c.parent_index)
                        if len(n.data.motif_name) > 0:
                            if n.data.motif_name[2] == "=":
                                ss = ssf.ss_id_to_secondary_structure(n.data.end_ss_id)
                                parent_end_ss_id = secondary_structure.assign_end_id(ss, ss.ends[1])
                    d = MotifTreeTopologyNodeData(c.name, c.end_id, parent_end_ss_id)
                    self.tree.add_data(d, c.parent_index)
                else:
                    if i > 0:
                        parent_end_ss_id = c.parent_end_id
                        n = self.tree.get_node(c.parent_index)
                        if len(n.data.motif_name) > 0:
                            if n.data.motif_name[2] == "=":
                                ss = ssf.ss_id_to_secondary_structure(n.data.end_ss_id)
                                parent_end_ss_id = secondary_structure.assign_end_id(ss, ss.ends[1])

                    ss = ssf.ss_id_to_secondary_structure(c.end_id)
                    conn = ss.motif_topology_from_end(ss.ends[0])
                    for j, c2 in enumerate(conn):
                        spl = c2[0].split("_")
                        name = spl[0][0]+spl[2][1]+"="+spl[0][1]+spl[2][0]
                        if j == 0:
                            d = MotifTreeTopologyNodeData(name, c2[0], c.parent_end_id)
                            self.tree.add_data(d, c.parent_index)
                        else:
                            d = MotifTreeTopologyNodeData(name, c2[0], c2[1])
                            self.tree.add_data(d)
            else:
                d = MotifTreeTopologyNodeData(c.name, c.end_id, c.parent_end_id)
                self.tree.add_data(d, c.parent_index)

            #print c
            #print "length",len(self.tree), self.tree.last_node.index


    def __repr__(self):
        s = "<(MotifTreeTopology: #nodes: %d\n" % (len(self.tree))
        for n in self.tree.nodes:

            s += "\t index: %d, name: %s, end_id: %s parent_end_id: %s, parent %s \n" %\
                 (n.index, n.data.motif_name, n.data.end_ss_id,
                  n.data.parent_end_ss_id, n.parent_index())
        s += ")"

        return s


class MotifTreeTopologyNodeData(object):
    def __init__(self, motif_name, end_ss_id, parent_end_ss_id):
        self.motif_name, self.end_ss_id = motif_name, end_ss_id
        self.parent_end_ss_id = parent_end_ss_id




