import os
import motif_scorer
import motif
import structure
import chain
import util
import pdb_parser
import x3dna
import residue
import basepair
import secondary_structure
import secondary_structure_factory
import settings


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

    def build_chains(self, residues):
        """
        takes all residues and puts into the correct order in chains checking
        for physical connection between O5' and P atoms between residues

        :param residues: residue objects that belong in this structure
        :type residues: List of Residue objects
        """

        chains = []
        # sort residues so check residues for connection quicker as the next on
        # in the array will be closest to it by number
        residues.sort(key=lambda x: x.num)

        while True:
            current = None
            # find next 5' end, all chains go from 5' to 3'
            for i, r in enumerate(residues):
                five_prime_end = 1
                for j, r2 in enumerate(residues):
                    if r.connected_to(r2) == -1:
                        five_prime_end = 0
                        break
                if five_prime_end:
                    current = r
                    break
            if not current:
                break
            residues.remove(current)
            current_chain_res = []
            # extend chain until 3' end
            while current is not None:
                current_chain_res.append(current)
                found = 0
                for r in residues:
                    if current.connected_to(r) == 1:
                        current = r
                        found = 1
                        break
                if found:
                    residues.remove(current)
                else:
                    # no more residues to add, make chain object
                    chains.append(chain.Chain(current_chain_res))
                    current = None

        return chains

    def get_structure(self, pdb_path):
        residues = pdb_parser.parse(pdb_path)
        chains = self.build_chains(residues)
        s = structure.Structure(chains)
        s.name = util.filename(pdb_path[:-4])
        return s

    def _setup_basepairs(self, path, name, structure):
        """
        gets x3dna data on basepairing information and then interwines it
        with the structural information stored in structure for simpler
        retreival of data
        """
        x3dna_parser = x3dna.X3dna()
        x_basepairs = x3dna_parser.get_basepairs(path)
        basepairs = []
        for xbp in x_basepairs:
            res1 = structure.get_residue(num=xbp.res1.num,
                                         chain_id=xbp.res1.chain_id,
                                         i_code=xbp.res1.i_code)

            res2 = structure.get_residue(num=xbp.res2.num,
                                         chain_id=xbp.res2.chain_id,
                                         i_code=xbp.res2.i_code)

            if res1 is None or res2 is None:
                raise ValueError("cannot find residues in basepair")

            bp = basepair.Basepair(res1, res2, xbp.r, xbp.bp_type)
            basepairs.append(bp)

        if os.path.isfile("ref_frames.dat"):
            os.remove("ref_frames.dat")

        if os.path.isfile(name + "_dssr.out"):
            os.remove(name + "_dssr.out")

        return basepairs

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
        best_i = -1
        for i, end in enumerate(m.ends):
            c1 = util.center(end.atoms)
            dist = util.distance(c1, c2)
            if dist < best:
                best = dist
                closest = end
                best_i = i

        updated_ends = [ closest ]
        for i, end in enumerate(m.ends):
            if end != closest:
                updated_ends.append(end)

        for i, end in enumerate(updated_ends):
            flip_res = 0
            for c in m.chains():
                if c.first() == end.res2:
                    flip_res = 1
                    break

            if flip_res:
                updated_ends[i].res1,updated_ends[i].res2 = \
                    updated_ends[i].res2,updated_ends[i].res1

        m.ends = updated_ends

    def _setup_secondary_structure(self, m):
        ss = secondary_structure_factory.factory.secondary_structure_from_motif(m)
        m.end_ids = ["" for x in m.ends]
        for i, end in enumerate(m.ends):
            res1 = ss.get_residue(uuid=end.res1.uuid)
            res2 = ss.get_residue(uuid=end.res2.uuid)
            ss_end = ss.get_bp(res1, res2)
            m.end_ids[i] = secondary_structure.assign_end_id(ss, ss_end)

        m.secondary_structure = ss

    def motif_from_file(self, path):
        filename = util.filename(path)
        # is a motif directory
        if os.path.isdir(path):
            structure = self.get_structure(path + "/" + filename + ".pdb")

        else:
            structure = self.get_structure(path)
            filename = filename[:-4]

        basepairs = self._setup_basepairs(path, filename, structure)
        ends = self._setup_basepair_ends(structure, basepairs)

        m           = motif.Motif()
        m.name      = filename
        m.path      = path
        m.structure = structure
        m.basepairs = basepairs
        m.ends      = ends
        m.score     = self.scorer.score(m)

        #try:
        self._setup_secondary_structure(m)
        #except:
        #    print "did not parse secondary_structure", m.name

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
        chains = self.build_chains(res)
        m.structure.chains = chains
        m.basepairs = bps
        m.ends = [bps[0], bps[-1]]
        self._setup_secondary_structure(m)
        return m

    def motif_from_res(self, res, bps):
        m = motif.Motif()
        chains = self.build_chains(res)
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


