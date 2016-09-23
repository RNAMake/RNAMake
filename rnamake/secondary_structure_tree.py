import tree

def tree_from_pose(p, start_m=None):
    sst = SecondaryStructureTree()
    for m in p.motifs:
        found = 0
        for m1 in p.motifs:
            if m1 == m:
                continue
            for e in m1.ends:
                if e == m.ends[0]:
                    found = 1
                    break
        if not found:
            start_m = m
            break

    if start_m is None:
        raise ValueError("cannot convert pose to tree there are no motifs with free ends")
    open_motifs = [[start_m, -1, -1]]
    seen_motifs = {}
    while len(open_motifs) > 0:
        current, parent_index, parent_end_index = open_motifs.pop(0)
        seen_motifs[current] = 1

        pos = sst.add_motif(current, parent_index, parent_end_index)
        for m in p.motifs:
            if m in seen_motifs:
                continue
            for i, end1 in enumerate(current.ends):
                if end1 == 0:
                    continue
                if end1 == m.ends[0]:
                    open_motifs.append([m, pos, i])
    return sst


class SecondaryStructureTree(object):
    def __init__(self):
        self.tree = tree.TreeStatic()

    def add_motif(self, m, parent_index=-1, parent_end_index=-1, m_end_index=0,
                  parent_end_name=None):
        parent = self.tree.last_node
        if parent_index != -1:
            parent = self.tree.get_node(parent_index)

        m_copy = m.copy()
        if parent is None:
            return self.tree.add_data(m_copy, len(m_copy.ends), -1, -1)

        avail_pos = self.tree.get_available_pos(parent, parent_end_index)
        for p in avail_pos:
            pos = self.tree.add_data(m_copy, len(m_copy.ends),
                                     parent.index, p)
            return pos

        return -1

