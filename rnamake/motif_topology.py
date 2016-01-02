
import motif_graph
import motif_tree
import motif_type
import resource_manager as rm

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
