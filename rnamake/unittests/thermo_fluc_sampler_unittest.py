import unittest
import random
import build
import math
import numpy as np
import numpy.linalg as lin
import rnamake.thermo_fluc_sampler
import rnamake.motif_tree_topology as motif_tree_topology
import rnamake.motif_state_ensemble_tree as motif_state_ensemble_tree
import rnamake.motif_tree as motif_tree
import rnamake.segmenter as segmenter
import rnamake.resource_manager as rm
import rnamake.pose_factory as pf
import rnamake.transformations as trans
import rnamake.util as util

class ThermoFlucSamplerUnittest(unittest.TestCase):

    def _fill_basepairs_in_ss(self, ss):
        pairs = ["AU", "UA", "GC", "CG"]
        for bp in ss.basepairs:
            if bp.res1.name == "N" and bp.res2.name == "N":
                p = random.choice(pairs)
                bp.res1.name = p[0]
                bp.res2.name = p[1]

    def test_creation(self):
        tfs = rnamake.thermo_fluc_sampler.ThermoFlucSampler()

    def test_sample(self):
        builder = build.BuildSecondaryStructure()
        ss = builder.build_helix(5)
        con = ss.motif_topology_from_end()
        mtt = motif_tree_topology.MotifTreeTopology(con)
        mt = motif_tree.motif_tree_from_topology(mtt)
        mset =  motif_state_ensemble_tree.MotifStateEnsembleTree(mt)

        tfs = rnamake.thermo_fluc_sampler.ThermoFlucSampler()
        tfs.setup(mset)
        tfs.sample(100)

    def test_relaxer(self):
        s = rnamake.segmenter.Segmenter()
        path = rnamake.settings.UNITTEST_PATH + "/resources/motifs/p4p6"
        p = pf.factory.pose_from_file(path)
        end1 = p.get_basepair(name='A111-A209')[0]
        end2 = p.get_basepair(name='A118-A203')[0]
        segments = s.apply(p, [end1, end2])
        pd = segments.remaining
        pd.name = "p4p6_frag"
        rm.manager.register_motif(pd)
        mt = motif_tree.MotifTree(sterics=0)
        mt.add_motif(pd)
        mt.add_motif(rm.manager.get_motif(name='HELIX.IDEAL.6'),
                     parent_end_name='A111-A209')
        mt.add_connection(0, 1, 'A118-A203')
        ss = mt.designable_secondary_structure()
        self._fill_basepairs_in_ss(ss)
        con = ss.motif_topology_from_end()
        mtt = motif_tree_topology.MotifTreeTopology(con)
        mt2 = motif_tree.motif_tree_from_topology(mtt, sterics=0)

        ni_1 = 0
        ei_1 = 0
        ni_2 = mt2.last_node().index
        ei_2 = 1
        mset =  motif_state_ensemble_tree.MotifStateEnsembleTree(mt2)

        relaxer = rnamake.thermo_fluc_sampler.ThermoFlucRelax()
        relaxer.run(mset, ni_1, ni_2, ei_1, ei_2)
        relaxer.write_pdbs()

    def test_update(self):
        builder = build.BuildSecondaryStructure()
        ss = builder.build_helix(5)
        con = ss.motif_topology_from_end()
        mtt = motif_tree_topology.MotifTreeTopology(con)
        mt = motif_tree.motif_tree_from_topology(mtt)
        mset =  motif_state_ensemble_tree.MotifStateEnsembleTree(mt)

        tfs = rnamake.thermo_fluc_sampler.ThermoFlucSampler()
        tfs.option('temperature', 1000.0)
        tfs.setup(mset)

        n = tfs.mst.get_node(3)
        for i in range(100):
            tfs.next()
            #print n.data.cur_state.end_states[0].d

    def _test_interia_rmsd(self):
        motif = rm.manager.get_motif(name='HELIX.IDEAL.3')
        Ix = 1.94917
        Iy = 12.3736
        Iz = 330.499

        r1 = motif.ends[0].r()
        r2 = motif.ends[1].r()
        r = util.unitarize(r1.T.dot(r2))
        e = trans.euler_from_matrix(r)

        x = math.atan2(r[2][1], r[2][2])
        y = math.atan2(-r[2][0], math.sqrt(r[2][1]*r[2][1] + r[2][2]*r[2][2]))
        z = math.atan2(r[1][0], r[0,0])

        g = e[2]
        b = e[1]
        a = e[0]


        t = motif.ends[0].d() - motif.ends[1].d()
        #t = [0, 0, 10]

        #g, b, a = 1.57, 0, 0

        k = (1 - (math.cos(g) * math.cos(a) * math.cos(b)) + (math.sin(a) * math.sin(g)) )
        l = (1 + (math.sin(g) * math.sin(a) * math.cos(b)) - (math.cos(a) * math.cos(g)) )
        m = (1 - math.cos(b))

        print k,l, m

        print Ix*k, Iy*l, Iz*m, lin.norm(t)
        rmsd = math.sqrt( 2 * ( (Ix * k) + (Iy * l) + (Iz * m) )  + lin.norm(t) ** 2 )

        #print  motif.ends[0].res1
        #print  motif.ends[1].res2

        dist_squared = 0
        count = 0
        for i, a in enumerate(motif.ends[0].res1.atoms):
            if a is None:
                continue
            dist_squared +=  util.distance(a.coords, motif.ends[1].res2.atoms[i].coords) ** 2
            count += 1
        for i, a in enumerate(motif.ends[0].res2.atoms):
            if a is None:
                continue
            dist_squared +=  util.distance(a.coords, motif.ends[1].res1.atoms[i].coords) ** 2
            count += 1
        dist_squared = dist_squared / count
        print rmsd
        print math.sqrt(dist_squared)



def main():
    unittest.main()

if __name__ == '__main__':
    main()