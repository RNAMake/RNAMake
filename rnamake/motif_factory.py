import os

import motif_scorer
import motif
import structure
import chain
import util
import x3dna
import residue
import basepair
import settings
import rna_structure
import motif_type
import pdb_parser
import exceptions
import residue_type
import transform
import bead

from rnamake.primitives.rna_structure import ends_from_basepairs, assign_end_id, end_id_to_seq_and_db

class MotifFactory(object):

    class _MotifElements(object):
        def __init__(self, structure, basepairs, ends):
            self.structure = structure
            self.basepairs = basepairs
            self.ends = ends

        @classmethod
        def copy(cls, me):
            s_new = structure.Structure.copy(me.structure)
            bps_new = []
            for bp in me.basepairs:
                bps_new.append(basepair.Basepair.copy(bp))
            ends_new = []
            for end in me.ends:
                i = me.basepairs.index(end)
                ends_new.append(bps_new[i])

            return cls(s_new, bps_new, ends_new)

    def __init__(self, rts=None):
        self._rts = rts
        if self._rts is None:
            self._rts = residue_type.ResidueTypeSet()

        self._clash_radius = settings.CLASH_RADIUS
        self._ref_motif = self.__setup_ref_motif()
        self._base_helix = self.__setup_base_motif()
        self._added_helix = motif.Motif.copy(self._base_helix)
        self._scorer = motif_scorer.MotifScorer()

    def __get_pdb_path(self, path):
        filename = util.filename(path)
        if os.path.isdir(path):
            return path + "/" + filename + ".pdb"

        else:
            return path

    def __setup_ref_motif(self):
        ref_motif_dir = settings.RESOURCES_PATH + "/start"
        pdb_path = self.__get_pdb_path(ref_motif_dir)
        s = structure.structure_from_pdb(pdb_path, self._rts)

        basepairs = rna_structure.basepairs_from_x3dna(pdb_path, s)
        ends      = ends_from_basepairs(s, basepairs)
        end_ids   = self.__setup_end_ids(s, basepairs, ends)

        name = "start"
        score = 0
        seq, dot_bracket = end_id_to_seq_and_db(end_ids[0])
        m = motif.Motif(s, basepairs, ends, end_ids, name, motif_type.HELIX, 0,
                        dot_bracket=dot_bracket)
        return m

    def __setup_base_motif(self):
        path = settings.RESOURCES_PATH + "base_helix"
        pdb_path = self.__get_pdb_path(path)
        s = structure.structure_from_pdb(pdb_path, self._rts)

        basepairs = rna_structure.basepairs_from_x3dna(pdb_path, s)
        ends      = ends_from_basepairs(s, basepairs)

        aligned_s    = self.__get_aligned_structure(s)
        aligned_ends = self.__get_aligned_ends(s, ends)
        end_ids      = self.__setup_end_ids(aligned_s, basepairs, aligned_ends)
        name = "Base"
        score = 0

        seq, dot_bracket = end_id_to_seq_and_db(end_ids[0])
        m = motif.Motif(s, basepairs, ends, end_ids, name, motif_type.HELIX, 0,
                        dot_bracket=dot_bracket)

        return m

    def __setup_end_ids(self, s, basepairs, ends):
        end_ids = []
        for end in ends:
            end_id = assign_end_id(s, basepairs, end)
            end_ids.append(end_id)
        return end_ids

    def __steric_clash(self, m, me):
        for r1 in m.iter_res():
            for b1 in r1.iter_beads():
                for r2 in me.structure.iter_res():
                    for b2 in r2.iter_beads():
                        if b1.btype == bead.BeadType.PHOS or \
                           b2.btype == bead.BeadType.PHOS:
                            continue
                        dist = b1.distance(b2)
                        if dist < self._clash_radius:
                            return 1

        return 0

    def __residue_overlap(self, m, me, ei1, ei2):
        end_res_1 = m.get_bp_res(m.get_end(ei1))
        end_res_2 = [
            me.structure.get_residue(uuid=me.ends[ei2].res1_uuid),
            me.structure.get_residue(uuid=me.ends[ei2].res2_uuid)
        ]

        for r1 in m.iter_res():
            if r1 in end_res_1:
                continue
            for r2 in me.structure.iter_res():
                if r2 in end_res_2:
                    continue
                dist = util.distance(r1.center(), r2.center())
                if dist < 3:
                    return 1

        return 0

    def __get_aligned_structure(self, s):
        closest = None
        best = 1000
        c2 = self._ref_motif.get_end(0).d
        for i, c in enumerate(s):
            c1 = c.first().center()
            dist = util.distance(c1, c2)
            if dist < best:
                best = dist
                closest = c
        updated_chains = [closest]
        for i, c in enumerate(s):
            if c != closest:
                updated_chains.append(c)

        return structure.Structure(updated_chains)

    def __get_aligned_ends(self, s, ends):
        if len(ends) == 0:
            return ends

        c2 = self._ref_motif.get_end(0).d
        closest = None
        best = 10000
        best_i = -1
        for i, end in enumerate(ends):
            c1 = end.d
            dist = util.distance(c1, c2)
            if dist < best:
                best = dist
                closest = end
                best_i = i
        updated_ends = [ closest ]
        for i, end in enumerate(ends):
            if end != closest:
                updated_ends.append(end)

        for i, end in enumerate(updated_ends):
            flip_res = 0
            for c in s:
                if c.first().uuid == end.res2_uuid:
                    flip_res = 1
                    break

            if flip_res:
                end.flip_res()

        return updated_ends

    def __align_motif_elements_to_frame(self, ref_bp, me, pos):
        me_copy = self._MotifElements.copy(me)

        r1 , r2 = ref_bp.r ,me_copy.ends[pos].r
        r = util.unitarize(r1.T.dot(r2))
        trans = -me_copy.ends[pos].d
        t = transform.Transform(r, trans)

        me_copy.structure.transform(t)
        for bp in me_copy.basepairs:
            bp.transform(t)

        bp_pos_diff = ref_bp.d - me_copy.ends[pos].d

        me_copy.structure.move(bp_pos_diff)
        for bp in me_copy.basepairs:
            bp.move(bp_pos_diff)

        # alignment is by center of basepair, it can be slightly improved by
        # aligning the c1' sugars
        res1_coord, res2_coord = me_copy.ends[pos].sugars
        ref_res1_coord, ref_res2_coord = ref_bp.sugars

        dist1 = util.distance(res1_coord, ref_res1_coord)
        dist2 = util.distance(res2_coord, ref_res1_coord)

        if dist1 < dist2:
            sugar_diff_1 = ref_res1_coord - res1_coord
            sugar_diff_2 = ref_res2_coord - res2_coord
        else:
            sugar_diff_1 = ref_res1_coord - res2_coord
            sugar_diff_2 = ref_res2_coord - res1_coord

        if dist1 < 5 or dist2 < 5:
            diff = (sugar_diff_1 + sugar_diff_2) / 2
            me_copy.structure.move(diff)
            for bp in me_copy.basepairs:
                bp.move(diff)

        for c in me_copy.structure:
            for r in c:
                found = 0
                for i, end in enumerate(me_copy.ends):
                    if (end.res1_uuid == r.uuid or end.res2_uuid == r.uuid) and \
                       i != pos:
                        found = 1
                        break
                if not found:
                    r.build_beads()

        return me_copy

    def __get_standardized_elements(self, me, ei):
        aligned_me = self.__align_motif_elements_to_frame(self._base_helix.get_end(1), me, ei)

        if self.__steric_clash(self._base_helix, aligned_me) or \
           self.__residue_overlap(self._base_helix, aligned_me, 1, ei):
            aligned_me.ends[ei].flip()
            aligned_me = self.__align_motif_elements_to_frame(self._base_helix.get_end(1), aligned_me, ei)

            if self.__steric_clash(self._base_helix, aligned_me) or \
               self.__residue_overlap(self._base_helix, aligned_me, 1, ei):
                return None

        fail = 0
        for i in range(len(aligned_me.ends)):
            if i == ei:
                continue

            m2_added = motif.get_aligned_motif(aligned_me.ends[i],
                                               self._added_helix.get_end(0),
                                               self._added_helix)

            if not self.__steric_clash(m2_added, aligned_me) and \
               not motif.clash_between_motifs(self._base_helix, m2_added):
                continue

            aligned_me.ends[i].flip()

            m2_added = motif.get_aligned_motif(aligned_me.ends[i],
                                               self._added_helix.get_end(0),
                                               self._added_helix)

            if self.__steric_clash(m2_added, aligned_me) or \
               motif.clash_between_motifs(self._base_helix, m2_added):
                fail = 1
                break

        if fail:
            return None
        else:
            return aligned_me

    def __filter_element_basepairs(self, element):
        basepairs = []

        for bp in element.basepairs:
            if bp.bp_type != "cW-W":
                basepairs.append(bp)
            if bp in element.ends:
                basepairs.append(bp)

        element.basepairs = basepairs

    def __motifs_from_elements(self, elements, mtype, m_name):
        score = self._scorer.score_elements(elements.structure, elements.basepairs)
        block_end_add = 0
        if mtype == motif_type.HAIRPIN:
            block_end_add=-1

        motifs = []
        for i in range(len(elements.ends)):
            aligned_me = self.__get_standardized_elements(elements, i)
            if aligned_me is None:
                continue

            aligned_s    = self.__get_aligned_structure(aligned_me.structure)
            aligned_ends = self.__get_aligned_ends(aligned_s, aligned_me.ends)

            end_bp = None
            for end in elements.ends:
                if end == aligned_ends[0]:
                    end_bp = end
                    break


            new_me = self._MotifElements(aligned_s, aligned_me.basepairs, aligned_ends)
            final_me = self.__align_motif_elements_to_frame(end_bp, new_me, 0)

            diff = util.distance(elements.structure._residues[0].center(),
                                final_me.structure._residues[0].center())

            if diff > 1:
                end_bp.flip()
                final_me = self.__align_motif_elements_to_frame(end_bp, new_me, 0)
                end_bp.flip()


            if mtype is not motif_type.HELIX and mtype is not motif_type.UNKNOWN:
                self.__filter_element_basepairs(final_me)

            end_ids      = self.__setup_end_ids(final_me.structure, final_me.basepairs,
                                                final_me.ends)

            if len(end_ids) > 0:
                seq, dot_bracket = end_id_to_seq_and_db(end_ids[0])
            else:
                dot_bracket = ""

            m = motif.Motif(final_me.structure, final_me.basepairs, final_me.ends, end_ids,
                            m_name, mtype, score, block_end_add=block_end_add,
                            dot_bracket=dot_bracket)

            motifs.append(m)

        return motifs

    def motifs_from_file(self, path, mtype=motif_type.UNKNOWN, include_protein=0):
        pdb_path = self.__get_pdb_path(path)
        s = structure.structure_from_pdb(pdb_path, self._rts)

        basepairs = rna_structure.basepairs_from_x3dna(pdb_path, s)
        ends      = ends_from_basepairs(s, basepairs)

        elements = self._MotifElements(s, basepairs, ends)
        name = util.filename(pdb_path)[:-4]

        return self.__motifs_from_elements(elements, mtype, name)


        """if include_protein:
            if os.path.isdir(path):
                p_residues = pdb_parser.parse(path + "/" + filename + ".pdb",
                                              protein=1, rna=0)
            else:
                p_residues = pdb_parser.parse(path, protein=1, rna=0)
            beads = []
            for r in p_residues:
                a = r.get_atom("CA")
                beads.append(residue.Bead(a.coords, residue.BeadType.BASE))
            m.protein_beads = beads

        return m"""

    def align_motif_to_common_frame(self, m, ei):
        return motif.get_aligned_motif(self._ref_motif.get_end(0), m.get_end(ei), m)

    def motifs_from_bps(self, bps, org_m, m_name, mtype):
        all_res = []
        for bp in bps:
            bp_res = org_m.get_bp_res(bp)
            for r in bp_res:
                if r not in all_res:
                    all_res.append(r)

        all_res = [residue.Residue.copy(r, build_beads=0) for r in all_res]

        chains = chain.connect_residues_into_chains(all_res)
        s = structure.Structure(chains)
        ends = ends_from_basepairs(s, bps)

        if len(bps) == 1:
            end_ids = self.__setup_end_ids(s, bps, ends)

            if len(end_ids) > 0:
                seq, dot_bracket = end_id_to_seq_and_db(end_ids[0])
            else:
                dot_bracket = ""

            m = motif.Motif(s, bps, ends, end_ids, m_name, mtype, 0,
                            block_end_add=-1, dot_bracket=dot_bracket)
            return m
        else:
            raise ValueError("not supported")


    def motifs_from_res(self, res, bps, org_m, m_name, mtype):
        all_res = []
        all_res.extend(res)
        for bp in bps:
            bp_res = org_m.get_bp_res(bp)
            for r in bp_res:
                if r not in all_res:
                    all_res.append(r)

        all_res = [ residue.Residue.copy(r, build_beads=0) for r in all_res]

        chains = chain.connect_residues_into_chains(all_res)
        s = structure.Structure(chains)
        ends = ends_from_basepairs(s, bps)
        elements = self._MotifElements(s, bps, ends)

        return self.__motifs_from_elements(elements, mtype, m_name)

    def __initial_aligned_ends(self, ends):
        if len(ends) == 0:
            return ends

        c2 = self._ref_motif.get_end(0).d
        closest = None
        best = 10000
        best_i = -1
        for i, end in enumerate(ends):
            c1 = end.d
            dist = util.distance(c1, c2)
            if dist < best:
                best = dist
                closest = end
                best_i = i
        updated_ends = [ closest ]
        for i, end in enumerate(ends):
            if end != closest:
                updated_ends.append(end)

        return updated_ends

    def rna_structure_from_element(self, s, bps, ends, name):

        elements = self._MotifElements(s, bps, ends)
        aligned_me = self.__align_motif_elements_to_frame(self._ref_motif.get_end(0), elements, 0)
        #aligned_me = elements

        #print s.get_chain(0).get_residue(0).name
        aligned_s = self.__get_aligned_structure(aligned_me.structure)
        aligned_ends = self.__get_aligned_ends(aligned_s, aligned_me.ends)
        #print aligned_s.get_chain(0).get_residue(0).name
        end_bp = None
        for end in elements.ends:
            if end == aligned_ends[0]:
                end_bp = end
                break

        new_me = self._MotifElements(aligned_s, aligned_me.basepairs, aligned_ends)
        final_me = self.__align_motif_elements_to_frame(end_bp, new_me, 0)
        final_me = new_me

        diff = util.distance(elements.structure._residues[0].center(),
                             final_me.structure._residues[0].center())

        if diff > 1:
            end_bp.flip()
            final_me = self.__align_motif_elements_to_frame(end_bp, new_me, 0)
            end_bp.flip()

        end_ids = self.__setup_end_ids(final_me.structure, final_me.basepairs,
                                       final_me.ends)

        if len(end_ids) > 0:
            seq, dot_bracket = end_id_to_seq_and_db(end_ids[0])
        else:
            dot_bracket = ""

        m = rna_structure.RNAStructure(final_me.structure, final_me.basepairs,
                                       final_me.ends, end_ids, name,
                                       dot_bracket=dot_bracket)

        return m



    def motifs_from_rstruc(self, rna_struc, mtype, name):
        chains = [ c for c in rna_struc.iter_chains() ]
        bps = [ bp for bp in rna_struc.iter_basepairs() ]
        ends = [ end for end in rna_struc.iter_ends() ]
        s = structure.Structure(chains)
        for r in s.iter_res():
            r.remove_beads()

        elements = self._MotifElements(s, bps, ends)
        return self.__motifs_from_elements(elements, mtype, name)

def ref_motif():
    path = settings.RESOURCES_PATH + "/start"
    m = factory.motif_from_file(path)
    return m



