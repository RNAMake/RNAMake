import util

class MTSCluster(object):
    def __init__(self, mts):
        self.state = mts.end_states[0]
        self.mtss = [mts]


class MotifCluster(object):
    def __init__(self, motif, state):
        self.state = state
        self.motifs = [motif]

#TODO dont just use first end_state
def cluster_mts(mtss, max_distance=1.5):
    clusters = [MTSCluster(mtss[0])]
    for i, mts in enumerate(mtss):
        if i == 0:
            continue
        found = 0
        for c in clusters:
            dist = c.state.diff(mts.end_states[0])
            if dist < max_distance:
                found = 1
                c.mtss.append(mts)
                break
        if not found:
            clusters.append(MTSCluster(mts))
    return clusters


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



