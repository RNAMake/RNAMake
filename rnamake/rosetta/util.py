from rnamake import motif_ensemble, motif_factory

import subprocess
import glob
import os

def convert_silent_file_to_motif_ensemble(path, keep=None):
    subprocess.call("extract_pdbs -in::file::silent %s 1000000" % (path), shell=True)

    pdbs = glob.glob("test_*.pdb")

    f = open(path)
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

    motifs = []
    kept_scores = []
    for j, pdb in enumerate(pdbs):
        m = motif_factory.factory.motif_from_file(pdb)
        os.remove(pdb)
        if keep is not None:
            res = []
            bps = []
            all_res = m.residues()
            for i in keep:
                res.append(all_res[i])
            for bp in m.basepairs:
                if bp.res1 in res and bp.res2 in res:
                    bps.append(bp)
            m_sub = motif_factory.factory.motif_from_res(res, bps)
            m = m_sub
        if len(m.ends) != 2:
            continue
        m_added = motif_factory.factory.can_align_motif_to_end(m, 0)
        if m_added is None:
            continue
        m_added = motif_factory.factory.align_motif_to_common_frame(m_added, 0)
        motifs.append(m_added)
        kept_scores.append(scores[j])

    me = motif_ensemble.MotifEnsemble()
    me.setup(motifs[0].end_ids[0], motifs, kept_scores)

    return me










