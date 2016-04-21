import os

import resource_manager as rm
import settings
import exceptions


class MotifClusters(object):
    def __init__(self):
        path = settings.RESOURCES_PATH + "/sim_list"
        if not os.path.isfile(path):
            raise exceptions.MotifClustersException("sim_list not in right place " +
                                                    "something very wrong ")
        self.sim_dict = {}

        f = open(path)
        lines = f.readlines()
        f.close()

        for l in lines:
            spl = l.split()
            key = spl[1]
            sims = spl[4:]
            self.sim_dict[key] = sims

    def get_clustered_motifs(self, m):
        if m.name not in self.sim_dict:
            raise exceptions.MotifClustersException("motif: " + m.name + " not"+
                                                    " not in clusters ")

        for m_name in self.sim_dict[m.name]:
            pass
