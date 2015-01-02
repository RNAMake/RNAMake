import util

class MTSCluster(object):
    def __init__(self, mts):
        self.state = mts.end_state
        self.mtss = [mts]

def cluster_mts(mtss, max_distance=1.5):
    clusters = [MTSCluster(mtss[0])]
    for i, mts in enumerate(mtss):
        if i == 0:
            continue
        found = 0
        for c in clusters:
            dist = c.state.diff(mts.end_state)
            if dist < max_distance:
                found = 1
                c.mtss.append(mts)
                break
        if not found:
            clusters.append(MTSCluster(mts))
    return clusters



