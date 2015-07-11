import motif_factory
import motif_type
import pose
import x3dna
import util

class PoseFactory(object):
    def __init__(self):
        pass

    def _copy_info_into_pose(self, m):
        p = pose.Pose()
        p.structure = m.structure
        p.ends      = m.ends
        p.end_ids   = m.end_ids
        p.basepairs = m.basepairs
        p.name      = m.name
        p.path      = m.path
        p.beads     = m.beads
        p.score     = m.score
        p.secondary_structure = m.secondary_structure

        return p

    #TODO implement finding tertiary contacts
    def _setup_motifs_from_x3dna(self, p, gu_are_helix, singlet_bp_seperation):
        x = x3dna.X3dna()
        x3dna_motifs = x.get_motifs(p.path)
        motifs = []
        for xm in x3dna_motifs:
            m = self._convert_x3dna_to_motif(xm, p)
            m.mtype = xm.mtype
            motifs.append(m)

        do_not_keep = {}
        pairs = []
        for i, m in enumerate(motifs):
            keep = 0

            if gu_are_helix and m.mtype == motif_type.TWOWAY:
                keep = self._remove_gu_basepair_motifs(m)
                if not keep:
                    do_not_keep[m] = 1

            if not singlet_bp_seperation and m.mtype == motif_type.TWOWAY:
                pair = self._merge_motifs_with_single_bp(m, motifs, i)
                if pair is not None:
                    do_not_keep[pair[0]] = 1
                    do_not_keep[pair[1]] = 1
                    pairs.append(pair)

        for m in motifs:
            if m not in do_not_keep:
                p.motif_dict[m.mtype].append(m)
                p.motif_dict[motif_type.ALL].append(m)

        for pair in pairs:
            res = []
            basepairs = []
            res.extend(pair[0].residues())
            basepairs.extend(pair[0].basepairs)
            for r in pair[1].residues():
                if r not in res:
                    res.append(r)
            for bp in pair[1].basepairs:
                if bp not in basepairs:
                    basepairs.append(bp)
            m = motif_factory.factory.motif_from_res(res, basepairs)
            m.mtype = motif_type.TWOWAY
            p.motif_dict[m.mtype].append(m)
            p.motif_dict[motif_type.ALL].append(m)




    def _remove_gu_basepair_motifs(self, m):
        res = []
        keep = 0
        for bp in m.basepairs:
            for r in bp.residues():
                if r not in res:
                    res.append(r)
            if not (util.wc_bp(bp) and bp.bp_type == "cW-W" or util.gu_bp(bp)):
                keep = 1

        if len(res) != len(m.residues()):
            keep = 1
        return keep

    def _merge_motifs_with_single_bp(self, m, motifs, i):
        for j, m2 in enumerate(motifs):
            if i >= j or m2.mtype != motif_type.TWOWAY:
                continue
            paired = 0
            for end in m.ends:
                if end in m2.ends:
                    return [m, m2]

            for c1 in m.chains():
                for c2 in m2.chains():
                    if abs(c1.first().num -c2.last().num) == 1:
                        return [m, m2]
                    if abs(c2.first().num -c1.last().num) == 1:
                        return [m, m2]

        return None


    def _convert_x3dna_to_motif(self, xm, p):
        res = []
        for xr in xm.residues:
            r = p.get_residue(num=xr.num, chain_id=xr.chain_id, i_code=xr.i_code)
            res.append(r)
        basepairs = []
        for r in res:
            bps = p.get_basepair(res1=r)
            for bp in bps:
                if bp.res1 in res and bp.res2 in res and bp not in basepairs:
                    basepairs.append(bp)

        m = motif_factory.factory.motif_from_res(res, basepairs)
        return m

    def pose_from_file(self, path, gu_are_helix=1, singlet_bp_seperation=0):
        base_motif = motif_factory.factory.motif_from_file(path)
        p = self._copy_info_into_pose(base_motif)
        self._setup_motifs_from_x3dna(p, gu_are_helix, singlet_bp_seperation)

        return p


    def pose_from_motif_tree(self):
        pass


factory = PoseFactory()