import util


class MotifCluster(object):
    def __init__(self, motif, state):
        self.state = state
        self.motifs = [motif]


def cluster_motifs(motifs, max_distance=1.5):
    clusters = [MotifCluster(motifs[0], motifs[0].ends[1].state())]
    for i, m in enumerate(motifs):
        if i == 0:
            continue
        found = 0
        for c in clusters:
            dist = c.state.diff(m.ends[1].state())
            if dist < max_distance:
                found = 1
                c.motifs.append(m)
                break
        if not found:
            clusters.append(MotifCluster(m, m.ends[1].state()))
    return clusters



