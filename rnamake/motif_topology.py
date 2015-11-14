
import motif_graph
import motif_tree
import motif_type
import resource_manager as rm

def graph_to_tree(mg, start=None):
        if start is None:
            start = mg.graph.oldest_node()

        seen = []
        open = [ start ]

        mt = motif_tree.MotifTree()
        mt.option('sterics', 0)
        ss = mg.secondary_structure()
        seen_nodes = {}
        index = 0

        for m in ss.motifs:
            print m
        exit()

        while len(open) > 0:
            current = open.pop(0)
            ss_m = ss.motif(current.data.id)

            if index == 0:
                free_end = -1
                for i, c in enumerate(current.connections):
                    if c is None:
                        free_end = i

                m = rm.manager.get_motif(name = ss_m.name,
                                         end_name=ss_m.ends[free_end].name(),
                                         end_id=ss_m.end_ids[free_end])

               # m = rm.manager.get_motif(name = ss_m.name)

                """try:
                    m = rm.manager.get_motif(end_name=ss_m.ends[free_end].name(),
                                             end_id=ss_m.end_ids[free_end])
                except:
                    m = rm.manager.get_motif(end_id=ss_m.end_ids[free_end])"""
                mt.add_motif(m)

            else:
                p = None
                p_index = None
                p_end_name = None
                c_end_index = -1
                c_end_name = None
                for c in current.connections:
                    if c is None:
                        continue
                    p_new = c.partner(current.index)
                    if p_new in seen_nodes:
                        p = p_new
                        p_index = seen_nodes[p_new]
                        p_end_index = c.end_index(p_new.index)
                        p_end_name = p.data.ends[p_end_index].name()
                        c_end_index = c.end_index(current.index)
                        c_end_name = current.data.ends[c_end_index].name()

                        break



                if ss_m.name[0:5] == "HELIX":
                    m = rm.manager.get_motif(name = ss_m.name)
                else:
                    m = rm.manager.get_motif(name = ss_m.name,
                                         end_name=ss_m.ends[free_end].name(),
                                         end_id=ss_m.end_ids[free_end])

                pos = mt.add_motif(m, parent_index=p_index, parent_end_name=p_end_name)
                if pos == -1:
                    raise ValueError("could not convert graph to tree")

            seen_nodes[current] = index
            index += 1
            for c in current.connections:
                if c is not None:
                    p = c.partner(current.index)
                    if p not in seen_nodes:
                        open.append(p)

        return mt
