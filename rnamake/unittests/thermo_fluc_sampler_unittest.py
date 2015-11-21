import unittest
import random
import build
import math
import numpy as np
import pandas as pd
import seaborn as sns
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
from rnamake import motif_topology, thermo_fluc_sampler, motif_graph

def center(points):
    length = points.shape[0]
    sum_x = np.sum(points[:, 0])
    sum_y = np.sum(points[:, 1])
    sum_z = np.sum(points[:, 2])

    return np.array([sum_x/length, sum_y/length, sum_z/length])

def rmsd(p1, p2):
    return np.linalg.norm(p1 - p2) / np.sqrt(len(p1))

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

    def get_moment_tensor(self, coords):
        sum = np.zeros([3,3])
        center = util.center_points(coords)
        for c in coords:
            diff = c - center
            c_T = diff[np.newaxis].T
            c = diff[np.newaxis]
            sum +=  c_T.dot(c)
        sum /= len(coords)
        eig, vec = lin.eig(sum)
        print eig
        print vec

    def get_moment(self, coords):
        center = util.center_points(coords)
        I = [0,0,0]
        for c in coords:
            I += (c - center) ** 2
        print I / len(coords)

    def get_moment_2(self, coords):
        center = util.center_points(coords)
        i_x = 0
        i_y = 0
        i_z = 0
        for c in coords:
            i_x += (c[1] - center[1] ) ** 2 + (c[2] - center[2]) ** 2
            i_y += (c[0] - center[0] ) ** 2 + (c[2] - center[2]) ** 2
            i_z += (c[0] - center[0] ) ** 2 + (c[1] - center[1]) ** 2

        return [i_x, i_y, i_z ]

    def test_interia_rmsd_3(self):
        motif = rm.manager.get_motif(name='HELIX.IDEAL')

        r1 = motif.ends[0].r()
        r2 = motif.ends[1].r()
        r = util.unitarize(r1.T.dot(r2))
        trans = lin.norm(motif.ends[0].base_d() - motif.ends[1].base_d())
        diff = motif.ends[0].base_d() - motif.ends[1].base_d()

        coords = []
        for i, a in enumerate( motif.ends[0].res1.atoms):
            if i < 12:
                continue
            coords.append(a.coords)
        for i, a in enumerate( motif.ends[0].res2.atoms):
            if i < 12:
                continue
            coords.append(a.coords)

        coords = np.array(coords)
        c = center(coords)
        #c = c / np.linalg.norm(c)

        dist_squared = 0
        count = 0
        for i, a in enumerate(motif.ends[0].res1.atoms):
            if a is None:
                continue
            if i < 12:
                continue
            dist_squared +=  util.distance(a.coords, motif.ends[1].res2.atoms[i].coords) ** 2
            count += 1
        for i, a in enumerate(motif.ends[0].res2.atoms):
            if a is None:
                continue
            if i < 12:
                continue
            dist_squared +=  util.distance(a.coords, motif.ends[1].res1.atoms[i].coords) ** 2
            count += 1
        dist_squared = dist_squared / count

        print rmsd(np.array([np.dot(c, r.T) + diff ]), np.array([c]))
        print math.sqrt(dist_squared)

    def test_interia_rmsd_2(self):
        motif = rm.manager.get_motif(name='HELIX.IDEAL.11')
        """Ix = 1.94917
        Iy = 12.3736
        Iz = 330.499"""

        #Ix= 6.20066117

        Ix = 9.606
        Iy = 1.94979491
        Iz = 0.02840045
        #Ix = 18.52865437
        #Iy = 3.048271
        #Iz = 21.525299703

        It = np.eye(3)
        It[0][0] = Ix
        It[1][1] = Iy
        It[2][2] = Iz
         #motif.ends[0].flip()
        r1 = motif.ends[0].r()
        r2 = motif.ends[1].r()
        r = util.unitarize(r1.T.dot(r2))
        e = trans.euler_from_matrix(r)
        g = e[2]
        b = e[1]
        a = e[0]

        t = lin.norm(motif.ends[0].base_d() - motif.ends[1].base_d())

        k = (1 - (math.cos(g) * math.cos(a) * math.cos(b)) + (math.sin(a) * math.sin(g)) )
        l = (1 + (math.sin(g) * math.sin(a) * math.cos(b)) - (math.cos(a) * math.cos(g)) )
        m = (1 - math.cos(b))

        #motif.ends[0].flip()
        r1 = motif.ends[0].r()
        r2 = motif.ends[1].r()
        r = util.unitarize(r1.T.dot(r2))
        t = lin.norm(motif.ends[0].base_d() - motif.ends[1].base_d())
        coords = []
        for i, a in enumerate( motif.ends[0].res1.atoms):
            if i < 12:
                continue
            coords.append(a.coords)
        for i, a in enumerate( motif.ends[0].res2.atoms):
            if i < 12:
                continue
            coords.append(a.coords)
        #self.get_moment_2(coords)
        #exit()
        #self.get_moment_tensor(coords)
        #exit()
        N = len(coords) - 10
        print math.sqrt(t**2), math.sqrt(2*(Ix + Iy + Iz)), math.sqrt(2*(Ix*k)), math.sqrt(2*(Iy*l)), math.sqrt(2*(Iz*m)), -math.sqrt(2*np.trace(r.dot(It)))
        #rmsd = math.sqrt(((t**2 + 2*(Ix +Iy + Iz) - 2*np.trace(r.dot(It))))/N)
        #rmsd = math.sqrt(t**2 + 2*(Ix +Iy + Iz) - 2*np.trace(r.dot(It)))
        #rmsd = math.sqrt(t**2 + 2*(Ix +Iy + Iz) - 2*np.trace(r.dot(It)))
        rmsd = math.sqrt(t**2 + 2*((Ix*k) + (Iy*l) + (Iz*m)))


        #print Ix*k, Iy*l, Iz*m, lin.norm(t)
        #rmsd = math.sqrt( 2 * ( (Ix * k) + (Iy * l) + (Iz * m) )  + lin.norm(t) ** 2 )

        #print  motif.ends[0].res1
        #print  motif.ends[1].res2

        dist_squared = 0
        count = 0
        for i, a in enumerate(motif.ends[0].res1.atoms):
            if a is None:
                continue
            if i < 12:
                continue
            dist_squared +=  util.distance(a.coords, motif.ends[1].res2.atoms[i].coords) ** 2
            count += 1
        for i, a in enumerate(motif.ends[0].res2.atoms):
            if a is None:
                continue
            if i < 12:
                continue
            dist_squared +=  util.distance(a.coords, motif.ends[1].res1.atoms[i].coords) ** 2
            count += 1
        dist_squared = dist_squared / count
        print rmsd
        #print math.sqrt(dist_squared) - math.sqrt(t**2)
        print math.sqrt(dist_squared)

    def test_interia_rmsd(self):
        motif = rm.manager.get_motif(name='HELIX.IDEAL')
        Ix = 18.52865437
        Iy = 3.048271
        Iz = 21.525299703

        r1 = motif.ends[0].r()
        r2 = motif.ends[1].r()
        r = util.unitarize(r1.T.dot(r2))
        e = trans.euler_from_matrix(r)

        g = e[2]
        b = e[1]
        a = e[0]

        #print r
        #print a,b,g


        t = lin.norm(motif.ends[0].base_d() - motif.ends[1].base_d())
        #t = [0, 0, 10]

        #g, b, a = 1.57, 0, 0

        k = (1 - (math.cos(g) * math.cos(a) * math.cos(b)) + (math.sin(a) * math.sin(g)) )
        l = (1 + (math.sin(g) * math.sin(a) * math.cos(b)) - (math.cos(a) * math.cos(g)) )
        m = (1 - math.cos(b))

        #print k,l, m
        #exit()

        print Ix*k, Iy*l, Iz*m, t
        rmsd = math.sqrt( 2 * ( (Ix * k) + (Iy * l) + (Iz * m) )  + t ** 2 )

        #print  motif.ends[0].res1
        #print  motif.ends[1].res2

        dist_squared = 0
        count = 0
        for i, a in enumerate(motif.ends[0].res1.atoms):
            if a is None:
                continue
            if i < 12:
                continue
            dist_squared +=  util.distance(a.coords, motif.ends[1].res2.atoms[i].coords) ** 2
            count += 1
        for i, a in enumerate(motif.ends[0].res2.atoms):
            if a is None:
                continue
            if i < 12:
                continue
            dist_squared +=  util.distance(a.coords, motif.ends[1].res1.atoms[i].coords) ** 2
            count += 1
        dist_squared = dist_squared / count
        print rmsd
        print math.sqrt(dist_squared)

    def test_folding(self):
        builder = build.BuildMotifGraph()
        #mg = builder.build(3)
        mg = motif_graph.MotifGraph()
        mg.add_motif(m_name="HELIX.IDEAL.12")
        mg.add_motif(m_name="TWOWAY.1GID.2")
        mg.add_motif(m_name="HELIX.IDEAL.12")
        mg.replace_ideal_helices()

        mt = motif_topology.graph_to_tree(mg, mg.last_node())
        mt.to_pdb("test.pdb", renumber=1)
        mset = motif_state_ensemble_tree.MotifStateEnsembleTree(mt=mt)

        tff = thermo_fluc_sampler.ThermoFlucFolding()
        tff.setup(mset, mt)
        #tff.movie = 1
        tff.run()

        df = pd.DataFrame(tff.contacts)
        sns.heatmap(df)
        sns.plt.show()







def main():
    unittest.main()

if __name__ == '__main__':
    main()