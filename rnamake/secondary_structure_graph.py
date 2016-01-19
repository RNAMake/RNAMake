import graph


def graph_from_pose(p, start_m=None):
    ssg = SecondaryStructureGraph()
    if start_m is None:
        start_m = p.motifs[0]
    seen_motifs = {}
    pos = -1

    free_end = None
    for end1 in start_m.ends:
        match = 0
        for m1 in p.motifs:
            if start_m == m1:
                continue
            for end2 in m1.ends:
                if end2 == end1:
                    match = 1
                    break
        if not match:
            free_end = end1
            break

    if free_end is None:
        raise ValueError("cannot convert pose to graph there is no free end to "
                         "start with")

    open_motifs = [[start_m, -1, -1, start_m.ends.index(free_end)]]
    while len(open_motifs) > 0:
        current, parent_index, parent_end_index, current_end_index = open_motifs.pop(0)
        seen_motifs[current] = 1

        pos = ssg.add_motif(current, parent_index, parent_end_index,
                            m_end_index=current_end_index)
        for m in p.motifs:
            if m in seen_motifs:
                continue
            for i, end1 in enumerate(current.ends):
                if i == current_end_index:
                    continue
                for j, end2 in enumerate(m.ends):
                    if end1 == end2:
                        open_motifs.append([m, pos, i, j])


    return ssg


class SecondaryStructureGraph(object):
    def __init__(self):
        self.graph = graph.GraphStatic()

    def __len__(self):
        return len(self.graph)

    def add_motif(self, m, parent_index=-1, parent_end_index=-1, m_end_index=0,
                  parent_end_name=None):

        parent = self.graph.last_node
        if parent_index != -1:
            parent = self.graph.get_node(parent_index)

        m_copy = m.copy()
        if parent is None:
            # do I need to create new res ids?
            pos = self.graph.add_data(m_copy, -1, -1, -1, len(m_copy.ends))
            return pos

        avail_pos = self.graph.get_availiable_pos(parent, parent_end_index)
        for p in avail_pos:
            pos = self.graph.add_data(m_copy, parent.index,
                                      p, m_end_index, len(m.ends))
            return pos



