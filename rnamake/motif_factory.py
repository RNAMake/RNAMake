import os

import motif_scorer
import motif
import structure
import chain
import util
import x3dna
import residue
import basepair
import secondary_structure
import secondary_structure_factory
import settings
import rna_structure
import motif_type
import pdb_parser
import exceptions


class MotifFactory(object):

    def __init__(self):
        try:
            path = settings.MOTIF_DIRS + "ref.motif"
            self.ref_motif = motif.file_to_motif(path)
            path = settings.MOTIF_DIRS + "base.motif"
            self.base_motif = motif.file_to_motif(path)
            self.base_motif.get_beads([self.base_motif.ends[1]])
            self.added_helix = self.base_motif.copy()
        except:
            pass
        self.clash_radius = settings.CLASH_RADIUS
        self.scorer = motif_scorer.MotifScorer()

    def _setup_basepair_ends(self, structure, basepairs):
        chain_ends = []
        for c in structure.chains:
            chain_ends.append(c.first())
            if len(c) > 1:
                chain_ends.append(c.last())

        ends = []
        for bp in basepairs:
            if bp.bp_type != "cW-W":
                continue
            if not (util.gu_bp(bp) or util.wc_bp(bp)):
                continue

            if bp.res1 in chain_ends and bp.res2 in chain_ends:
                ends.append(bp)

        return ends

    def _steric_clash(self, m1, m2):
        for c1 in m1.beads:
            for c2 in m2.beads:
                if c1.btype == residue.BeadType.PHOS or \
                   c2.btype == residue.BeadType.PHOS:
                    continue
                dist = util.distance(c1.center, c2.center)
                if dist < self.clash_radius:
                    return 1
        return 0

    def _align_chains(self, m):
        chains = m.chains()
        closest = None
        best = 1000
        best_i = -1
        c2 = util.center(self.ref_motif.ends[0].atoms)
        for i, c in enumerate(chains):
            c1 = util.center(c.first().atoms)
            dist = util.distance(c1, c2)
            if dist < best:
                best = dist
                closest = c
                best_i = i

            updated_chains = [closest]
            for i, c in enumerate(chains):
                if c != closest:
                    updated_chains.append(c)

            m.structure.chains = updated_chains

    def _align_ends(self, m):
        c2 = util.center(self.ref_motif.ends[0].atoms)
        closest = None
        best = 10000
        for i, end in enumerate(m.ends):
            c1 = util.center(end.atoms)
            dist = util.distance(c1, c2)
            if dist < best:
                best = dist
                closest = end

        updated_ends = [ closest ]
        for i, end in enumerate(m.ends):
            if end != closest:
                updated_ends.append(end)

        for i, end in enumerate(updated_ends):
            flip_res = 0
            for c in m.chains():
                if c.first().uuid == end.res2.uuid:
                    flip_res = 1
                    break
                if c.last().uuid == end.res1.uuid:
                    flip_res = 1
                    break

            #print end.name(), flip_res

            if flip_res:
                updated_ends[i].res1, updated_ends[i].res2 = \
                    updated_ends[i].res2, updated_ends[i].res1

        m.ends = updated_ends

    def _setup_secondary_structure(self, m):
        ss = secondary_structure_factory.factory.secondary_structure_from_motif(m)
        ss.end_ids = ["" for x in m.ends]
        m.end_ids = ["" for x in m.ends]
        #print m.name, len(m.ends)
        for i, end in enumerate(m.ends):
            res1 = ss.get_residue(uuid=end.res1.uuid)
            res2 = ss.get_residue(uuid=end.res2.uuid)
            ss_end = ss.get_basepair(res1, res2)
            m.end_ids[i] = secondary_structure.assign_end_id_new(ss, ss_end)
            ss.end_ids[i] = m.end_ids[i]

        m.secondary_structure = ss

    def motif_from_file(self, path, include_protein=0):
        filename = util.filename(path)
        # is a motif directory
        if os.path.isdir(path):
            s = structure.structure_from_pdb(path + "/" + filename + ".pdb")

        else:
            s = structure.structure_from_pdb(path)
            filename = filename[:-4]

        basepairs = rna_structure.basepairs_from_x3dna(path, filename, s)
        ends      = rna_structure.ends_from_basepairs(s, basepairs)

        r_struct = rna_structure.RNAStructure(s, basepairs, ends, filename,
                                              path, motif_type.UNKNOWN)

        #if len(r_struct.residues()) == 0:
        #    raise exceptions.MotifFactoryException(
        #   )

        m           = motif.Motif(r_struct)
        m.score     = self.scorer.score(m)

        #try:
        self._setup_secondary_structure(m)
        #except:
        #    print "did not parse secondary_structure", m.name

        if include_protein:
            p_residues = pdb_parser.parse(path, protein=1, rna=0)
            beads = []
            for r in p_residues:
                a = r.get_atom("CA")
                beads.append(residue.Bead(a.coords, residue.BeadType.BASE))
            m.protein_beads = beads

        return m

    def can_align_motif_to_end(self, m, ei):
        m_added = motif.get_aligned_motif(self.base_motif.ends[1], m.ends[ei], m)

        if self._steric_clash(self.base_motif, m_added):
            m.ends[ei].flip()
            m_added = motif.get_aligned_motif(self.base_motif.ends[1], m.ends[ei], m)

            if self._steric_clash(self.base_motif, m_added):
                return None

        fail = 0
        for i in range(len(m_added.ends)):
            if i == ei:
                continue
            m2_added = motif.get_aligned_motif(m_added.ends[i],
                                               self.added_helix.ends[0],
                                               self.added_helix)

            if not self._steric_clash(m_added, m2_added) and \
               not self._steric_clash(self.base_motif, m2_added):
                continue


            m_added.ends[i].flip()
            m.ends[i].flip()

            m2_added = motif.get_aligned_motif(m_added.ends[i],
                                               self.added_helix.ends[0],
                                               self.added_helix)

            if self._steric_clash(m_added, m2_added) or \
               self._steric_clash(self.base_motif, m2_added):
                fail = 1
                break

        if fail:
            return None
        else:
            return m_added

    def align_motif_to_common_frame(self, m, ei):
        m_added = motif.get_aligned_motif(self.ref_motif.ends[0], m.ends[ei], m)
        self.standardize_motif(m_added)
        return m_added

    def standardize_motif(self, m):
        self._align_chains(m)
        self._align_ends(m)
        try:
            self._setup_secondary_structure(m)
        except:
            print m.name
            return None

    def motif_from_bps(self, bps):
        m = motif.Motif()
        res = []
        for bp in bps:
            res.extend(bp.residues())
        chains = chain.connect_residues_into_chains(res)
        m.structure.chains = chains
        m.basepairs = bps
        m.ends = self._setup_basepair_ends(m.structure, m.basepairs)
        self._setup_secondary_structure(m)
        return m

    def motif_from_res(self, res, bps):
        m = motif.Motif()
        chains = chain.connect_residues_into_chains(res)
        m.structure.chains = chains
        m.basepairs = bps
        ends = self._setup_basepair_ends(m.structure, bps)
        m.ends = ends
        try:
            self._setup_secondary_structure(m)
        except:
            pass
        return m

    def motif_from_chains(self, chains, bps):
        m = motif.Motif()
        m.structure.chains = chains
        m.basepairs = bps
        ends = self._setup_basepair_ends(m.structure, bps)
        m.ends = ends
        try:
            self._setup_secondary_structure(m)
        except:
            pass
        return m



def ref_motif():
    path = settings.RESOURCES_PATH + "/start"
    m = factory.motif_from_file(path)
    return m

factory = MotifFactory()


