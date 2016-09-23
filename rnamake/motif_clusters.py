import os

import resource_manager as rm
import settings
import exceptions


class MotifClusters(object):
    def __init__(self):
        path = settings.RESOURCES_PATH + "/sim_list_new"
        if not os.path.isfile(path):
            raise exceptions.MotifClustersException("sim_list not in right place " +
                                                    "something very wrong ")
        self.sim_dict = {}

        f = open(path)
        lines = f.readlines()
        f.close()

        for l in lines:
            spl = l.split("|")
            key = spl[0].rstrip()
            sims = spl[1:-1]
            self.sim_dict[key] = sims

    def get_clustered_motifs(self, m):
        key = m.name + "," + m.ends[0].name()
        if key not in self.sim_dict:
            raise exceptions.MotifClustersException("motif: " + m.name + " not"+
                                                    " not in clusters ")
        motifs = []
        for motif_info in self.sim_dict[key]:
            spl = motif_info.split(",")
            new_m = rm.manager.get_motif(name=spl[0].lstrip(), end_name=spl[1].rstrip())
            motifs.append(new_m)
        return motifs
