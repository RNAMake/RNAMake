import subprocess
import glob
import os
from rnamake import motif_factory
from rnamake import motif_ensemble
from rnamake.submit import qsub_job

class RNADenovo(object):
    def __init__(self, nstruct=2):
        self.nstruct = nstruct

    def run(self, ss):
        cutpoints = []
        seq = ""
        db = ""
        length = 0
        for i, c in enumerate(ss.chains):
            seq += c.sequence()
            db  += c.dot_bracket()
            if i < len(ss.chains)-1:
                length += len(c.sequence())
                cutpoints.append(str(length))

        for r in ss.residues():
            print r.num, r.name,
        print

        subprocess.call('rna_denovo_setup.py -sequence \"%s\" -secstruct \"%s\" -cutpoint_open %s -out_script default.run -tag default -nstruct %d' % (seq, db, " ".join(cutpoints), self.nstruct), shell=True)

        subprocess.call("source default.run", shell=True)

        subprocess.call("extract_lowscore_decoys.py default.out " + str(self.nstruct),
                        shell=True)


    def generate_script(self, ss, name="test.sh", walltime="1:00:00"):
        cutpoints = []
        seq = ""
        db = ""
        length = 0
        for i, c in enumerate(ss.chains):
            seq += c.sequence()
            db  += c.dot_bracket()
            if i < len(ss.chains)-1:
                length += len(c.sequence())
                cutpoints.append(str(length))

        job_str = 'rna_denovo_setup.py -sequence \"%s\" -secstruct \"%s\" -cutpoint_open %s -out_script default.run -tag default -nstruct %d\n' % (seq, db, " ".join(cutpoints), self.nstruct)
        job_str += "source default.run\n"
        job_str += "extract_lowscore_decoys.py default.out " + str(self.nstruct)

        return qsub_job.QSUBJob(name, job_str, walltime)



    def process(self, ss, start_bp, name="default"):
        pdbs = glob.glob("default.out.*.pdb")

        f = open("default.out")
        lines = f.readlines()
        f.close()

        scores = []
        for l in lines:
            if l[0:5] != "SCORE":
                continue
            spl = l.split()
            try:
                scores.append(float(spl[1]))
            except:
                pass

        ss_res = ss.residues()
        motifs = []
        for pdb in pdbs:
            m = motif_factory.factory.motif_from_file(pdb)
            for i, r in enumerate(m.residues()):
                r.num = ss_res[i].num
                r.chain_id = ss_res[i].chain_id

            ei = m.get_end_index(start_bp.name())
            m_added = motif_factory.factory.can_align_motif_to_end(m, ei)
            if m_added is None:
                continue
            m_added = motif_factory.factory.align_motif_to_common_frame(m_added, ei)
            motifs.append(m_added)

        me = motif_ensemble.MotifEnsemble()
        print scores
        me.setup(motifs[0].end_ids[0], motifs, scores)
        s = me.to_str()
        f = open(name+".me", "w")
        f.write(s)
        f.close()

    def clean_up(self):
        os.remove("default.out")
        pdbs = glob.glob("default.out.*.pdb")
        for p in pdbs:
            os.remove(p)


