import collections

import motif_graph
import motif_tree
import motif_type
import motif
import resource_manager as rm
import exceptions


class GraphtoTree(object):

    def __init__(self):
        self.mt = motif_tree.MotifTree()

    class _GraphtoTreeNode(object):
        def __init__(self, parent, parent_end_index, node):
            self.parent = parent
            self.parent_end_index = parent_end_index
            self.node = node
            self.motif = node.data.copy()

    def _get_start_node(self, mg, start, start_end_index):
        if start is None:
            not_aligned = mg.get_not_aligned_nodes()
            if len(not_aligned) == 0 :
                raise RuntimeError("cannot convert graph to tree no starting point")
            start = not_aligned[0]
            start_n = self._GraphtoTreeNode(None, None, start)
        else:
            start_n = self._GraphtoTreeNode(None, None, start)
            if start_end_index is not None:
                new_motif = self._reorient_motif(start_n.motif, start_end_index)
                start_n.motif = new_motif

        return start_n

    def _get_new_nodes(self, current):
        new_nodes = []
        for c in current.node.connections:
            if c is None:
                continue
            partner = c.partner(current.node.index)
            # check to see if motif is already in tree
            try:
                n = self.mt.get_node_by_id(partner.data.id)
            except:
                parent_end_index = c.end_index(current.node.index)
                node_end_index = c.end_index(partner.index)
                new_n = self._GraphtoTreeNode(current.node, parent_end_index,
                                              partner)
                new_n.motif = self._reorient_motif(new_n.motif, node_end_index)
                new_n.parent_end_index = self._get_new_parent_end_index(current.node, c)
                new_nodes.append(new_n)

        return new_nodes

    def _get_new_parent_end_index(self, parent, c):

        if parent.data.mtype != motif_type.HELIX:
            parent_end_index = c.end_index(parent.index)
            parent_end_name = parent.data.ends[parent_end_index].name()
            tree_parent = self.mt.get_node_by_id(parent.data.id)

            for i, end in enumerate(tree_parent.data.ends):
                if end.name() == parent_end_name:
                    return i

            raise RuntimeError("did not find original end something went really wrong")

        # helices always go end 0 to 1
        else:
            return 1

    def _reorient_motif(self, m, end_index):
        # nothing needs to be done
        if end_index == 0:
            return m

        if m.mtype != motif_type.HELIX:
            try:
                new_m = rm.manager.get_motif(name=m.name,
                                             end_name=m.ends[end_index].name())

            except exceptions.ResourceManagerException as e:
                raise exceptions.MotifTopologyException(
                    "cannot convert graph to tree because: " + e.message)

            new_m.copy_uuids_from_motif(m)
            return new_m

        else:
            if m.name[:5] == "HELIX":
                new_m = rm.manager.get_motif(name=m.name)
                new_m.copy_uuids_from_motif(m)
                return new_m

            else:
                spl = m.name.split("=")
                new_name = spl[1] + "=" + spl[0]
                new_m = rm.manager.get_motif(name=new_name)
                new_m.id = m.id
                return new_m

    def convert(self, mg, start=None, start_end_index=None, last_node=None):
        self.mt = motif_tree.MotifTree()
        self.mt.option('sterics', 0)

        start_n = self._get_start_node(mg, start, start_end_index)

        open_nodes = [ start_n ]
        last_node_to_add = None

        while len(open_nodes) > 0:
            current = open_nodes.pop(0)

            if last_node == current.node:
                last_node_to_add = current
                continue

            if current.parent is None:
                self.mt.add_motif(current.motif)
            else:
                new_parent_index = self.mt.get_node_by_id(
                                    current.parent.data.id).index
                self.mt.add_motif(current.motif, parent_index=new_parent_index,
                                  parent_end_index = current.parent_end_index)
            new_nodes = self._get_new_nodes(current)

            # check for duplicates
            for n1 in new_nodes:
                found = 0
                for n2 in open_nodes:
                    if n2.node.index == n1.node.index:
                        found = 1
                        break
                    if last_node_to_add is not None:
                        if last_node_to_add.node.index == n1.index:
                            found = 1
                            break

                if not found:
                    open_nodes.append(n1)


        if last_node_to_add is not None:
            new_parent_index = self.mt.get_node_by_id(
                last_node_to_add.parent.data.id).index
            self.mt.add_motif(last_node_to_add.motif,
                              parent_index=new_parent_index,
                              parent_end_index=last_node_to_add.parent_end_index)

        # add missing connections
        for n in mg:
            for c in n.connections:
                if c is None:
                    continue
                n_id = n.data.id
                partner = c.partner(n.index)
                pn_id = partner.data.id

                try:
                    mt_n1 = self.mt.get_node_by_id(n_id)
                    mt_n2 = self.mt.get_node_by_id(pn_id)
                except exceptions.MotifTreeException:
                    continue

                if mt_n1.parent == mt_n2:
                    continue
                if mt_n2.parent == mt_n1:
                    continue

                n1_ei = c.end_index(n.index)
                if n.data.mtype == motif_type.HELIX:
                    n1_ei = 1

                n1_end_name = n.data.ends[n1_ei].name()

                if self.mt.connections.in_connection(mt_n1.index, n1_end_name):
                    continue


                n2_ei = c.end_index(partner.index)
                if partner.data.mtype == motif_type.HELIX:
                    n2_ei = 1

                n2_end_name = partner.data.ends[n2_ei].name()
                self.mt.add_connection(mt_n1.index, mt_n2.index,
                                      n1_end_name, n2_end_name)



        self.mt.option('sterics', 1)
        return self.mt



def graph_to_tree(mg, start=None, last_end=None):
        if start is None:
            start = mg.graph.oldest_node()

        seen = []
        open = [ start ]

        mt = motif_tree.MotifTree()
        mt.option('sterics', 0)
        ss = mg.secondary_structure()
        seen_nodes = {}
        seen_connections = {}
        index = 0
        pos = 0

        while len(open) > 0:
            current = open.pop(0)
            ss_m = ss.motif(current.data.id)

            if index == 0:
                free_end = -1
                for i, c in enumerate(current.connections):
                    if c is None:
                        free_end = i

                name = ss_m.name
                if ss_m.name[2] == "=":
                   m =  rm.manager.get_motif(name = name)
                else:
                    m = rm.manager.get_motif(name = name,
                                            end_name=ss_m.ends[free_end].name())
                mt.add_motif(m)

            else:
                highest = -1
                p = None
                p_index = None
                p_end_index = None
                p_end_name = None
                c_end_index = -1
                c_end_name = None
                for c in current.connections:
                    if c is None:
                        continue
                    p_new = c.partner(current.index)
                    if p_new in seen_nodes:
                        p_index = seen_nodes[p_new]
                        if highest < p_index:
                            highest = p_index
                        else:
                            continue
                        p = p_new
                        p_end_index = c.end_index(p_new.index)
                        p_end_name = p.data.ends[p_end_index].name()
                        c_end_index = c.end_index(current.index)
                        c_end_name = current.data.ends[c_end_index].name()
                        if p.data.name[2] == "=":
                            act_parent = mt.get_node(p_index)
                            p_end_index = 1
                            p_end_name = act_parent.data.ends[1].name()

                if last_end is not None:
                    if last_end.name() == p_end_name and len(seen_nodes) != len(mg)-1:
                        open.append(current)
                        continue


                if ss_m.name[2] == "=":
                    name = ss_m.name
                    if c_end_index == 1:
                        spl = name.split("=")
                        name = spl[1] + "=" + spl[0]
                    m = rm.manager.get_motif(name = name)

                else:
                    m = rm.manager.get_motif(name = ss_m.name,
                                             end_name=c_end_name)

                #print m.name, p.data.name, p_index, p_end_index, p_end_name
                seen_connections[str(p_index) + " " + str(index)] = 1
                pos = mt.add_motif(m, parent_index=p_index, parent_end_name=p_end_name)
                if pos == -1:
                    raise ValueError("could not convert graph to tree")

            seen_nodes[current] = index
            index += 1
            for c in current.connections:
                if c is not None:
                    p = c.partner(current.index)
                    if p not in seen_nodes and p not in open:
                        open.append(p)

        for n in mg.graph:
            for c in n.connections:
                if c is None:
                    continue
                key1 = str(seen_nodes[c.node_1]) + " " + str(seen_nodes[c.node_2])
                key2 = str(seen_nodes[c.node_2]) + " " + str(seen_nodes[c.node_1])
                if key1 not in seen_connections and key2 not in seen_connections:
                    mt.add_connection(seen_nodes[c.node_1], seen_nodes[c.node_2])
                    seen_connections[key1] = 1

        return mt
