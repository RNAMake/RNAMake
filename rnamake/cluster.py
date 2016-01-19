import util
import motif_ensemble
import math


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

def cluster_motifs_2(motifs, max_size=100):
    i = 0
    clusters = []
    for m in motifs:
        clusters.append(MotifCluster(m, m.ends[1].state()))
        i += 1
        if i == max_size:
            break

    if i < max_size:
        return clusters

    i = -1
    for m in motifs:
        i += 1
        if i < max_size:
            continue
        best_score = 10000
        best_c = None
        for c in clusters:
            dist = c.state.diff(m.ends[1].state())
            if dist < best_score:
                best_score = dist
                best_c = c
        best_c.motifs.append(m)
    return clusters

def cluster_motifs_3(motifs, start_dist=0.65,ideal_size=100):

    clusters = []
    if len(motifs) < ideal_size:
        for m in motifs:
            clusters.append(MotifCluster(m, m.ends[1].state()))
        return clusters

    start_dist = 0.65
    cluster_num = 0
    count = 1
    while abs(cluster_num - ideal_size) > 5:
        clusters = [MotifCluster(motifs[0], motifs[0].ends[1].state())]
        for i, m in enumerate(motifs):
            if i == 0:
                continue
            found = 0
            for c in clusters:
                dist = c.state.diff(m.ends[1].state())
                if dist < start_dist:
                    found = 1
                    c.motifs.append(m)
                    break
            if not found:
                clusters.append(MotifCluster(m, m.ends[1].state()))

        if count > 50:
            break

        cluster_num = len(clusters)
        if cluster_num - ideal_size < 0:
            start_dist -= 0.20 / count
        else:
            start_dist += 0.20 / count

        print cluster_num

        count += 1
    return clusters

def clusters_to_motif_ensemble(clusters, name="motif"):
    motifs = []
    energies = []

    kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
    kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius

    all_pop = 0
    for i, c_motifs in enumerate(clusters):
        all_pop += len(c_motifs.motifs)

    for i, c_motifs in enumerate(clusters):
        m = c_motifs.motifs[0]
        m.name = name + "." + str(i)
        pop =  float(len(c_motifs.motifs)) / float(all_pop)
        energy = -kBT*math.log(pop)
        motifs.append(m)
        energies.append(energy)

    me = motif_ensemble.MotifEnsemble()
    me.setup(motifs[0].end_ids[0], motifs, energies)
    return me


