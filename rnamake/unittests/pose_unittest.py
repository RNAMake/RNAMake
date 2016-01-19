import unittest
import build
import os
import rnamake.pose
import rnamake.settings as settings
import rnamake.motif_type as motif_type
import rnamake.motif_tree as motif_tree
import random
import rnamake.eternabot.sequence_designer as sequence_designer
import rnamake.pose_factory as pf
import rnamake.motif_factory as mf
from rnamake import motif_state_ensemble_tree
from rnamake import thermo_fluc_sampler
from rnamake import motif_outputer
from rnamake import motif_topology
from rnamake import resource_manager as rm
from rnamake import motif_ensemble
from rnamake.rosetta import rna_denovo

class PoseUnittest(unittest.TestCase):

    def test_creation(self):
        p = pf.factory.pose_from_file("resources/motifs/p4p6")

    def _test_designable_secondary_structure(self):
        builder = build.BuildMotifTree()
        mt = builder.build()
        p = mt.to_pose()
        ss = p.designable_secondary_structure()
        designer = sequence_designer.SequenceDesigner()
        results = designer.design(ss.dot_bracket(), ss.sequence())
        #print results[0]['end'][0]

    def test_motifs(self):
        p = pf.factory.pose_from_file("resources/motifs/p4p6")
        twoways = p.motifs(motif_type.TWOWAY)
        if len(twoways) != 6:
            self.fail("did not properly get all two way junctions")

    def _test_load_as_mg(self):
        p = pf.factory.pose_from_file("/Users/josephyesselman/projects/RNAMake.projects/claw/tRNA.pdb")

        all_motifs = p.motifs(motif_type.ALL)
        not_helix= []
        helix = []
        for m in all_motifs:
            if m.mtype == motif_type.HELIX:
                helix.append(m)
            else:
                not_helix.append(m)

        basepair_steps = []
        for m in helix:
            c = m.chains()[0]
            for i in range(1,len(c.residues)):
                sub_res = [c.residues[i-1], c.residues[i] ]
                sub_bps = []
                for r in sub_res:
                    bp = m.get_basepair(res1=r)[0]
                    for bp_r in bp.residues():
                        if bp_r not in sub_res:
                            sub_res.append(bp_r)
                    if bp not in sub_bps:
                         sub_bps.append(bp)

                bp_step = mf.factory.motif_from_res(sub_res, sub_bps)
                bp_step.mtype = motif_type.HELIX
                basepair_steps.append(bp_step)

        start = basepair_steps.pop(0)
        mt = motif_tree.MotifTree()
        all_motifs = basepair_steps + not_helix
        mt.add_motif(start)
        while len(all_motifs) > 0:
            leafs = mt.leafs_and_ends()
            for l in leafs:
                bp = l[0].data.ends[l[1]]
                found = 0
                next = None
                #print l[0].index, bp
                for m in all_motifs:
                    old_bps = m.get_basepair(name=bp.name())
                    if len(old_bps) == 0:
                        continue
                    pos = m.ends.index(old_bps[0])
                    next = m
                    if pos != 0:
                        print "failed"
                        exit()
                mt.write_pdbs("test")
                mt.add_motif(next, parent_index=l[0].index, parent_end_index=l[1])
                all_motifs.remove(next)
        mset = motif_state_ensemble_tree.MotifStateEnsembleTree(mt)
        sampler = thermo_fluc_sampler.ThermoFlucSampler()
        sampler.option('temperature', 1000.0)
        sampler.setup(mset)

        f = open("movie.pdb", "w")
        f.write("MODEL 1\n")
        f.write(sampler.to_pdb_str())
        f.write("ENDMDL\n")
        for i in range(100):
            if sampler.next() == 0:
                continue
            f.write("MODEL "+str(i+1)+"\n")
            f.write(sampler.to_pdb_str())
            f.write("ENDMDL\n")
        f.close()

    def test_load_as_mg2(self):
         #p = pf.factory.pose_from_file_new("/Users/josephyesselman/projects/RNAMake.projects/claw/tRNA.pdb")
         #p = pf.factory.pose_from_file_new("/Users/josephyesselman/projects/REDESIGN/examples/example_2/p4p6/p4p6.pdb")
        p = pf.factory.pose_from_file_new("/Users/josephyesselman/projects/DasLabBot/tests/runs/v1/1LNG/1LNG.pdb")

        #rna_denovo_runner = rna_denovo.RNADenovo(nstruct=5)
        for n in p.mgraph.graph:
            rm.manager.add_motif(motif=n.data)


        mt = motif_topology.graph_to_tree(p.mgraph)
        return

        for n in mt:
            if n.data.mtype != motif_type.HELIX:
                if n.index != 2:
                    continue
                me_name = "node." + str(n.index)
                if os.path.isfile(me_name):
                    continue
                print n.index, n.data.ends[0], n.data.sequence(), n.data.dot_bracket()
                #rna_denovo_runner.run(n.data.secondary_structure)
                #rna_denovo_runner.process(n.data.secondary_structure, n.data.ends[0],
                #                          me_name)
                #rna_denovo_runner.clean_up()
                #exit()

        mset = motif_state_ensemble_tree.MotifStateEnsembleTree(mt)
        sampler = thermo_fluc_sampler.ThermoFlucSampler()
        sampler.option('temperature', 1000.0)
        sampler.setup(mset)
        mt2 = sampler.mst.to_motif_tree()
        mt2.merger.to_pdb("test.pdb",renumber=1)





def main():
    unittest.main()

if __name__ == '__main__':
    main()
