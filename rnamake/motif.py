import itertools
import numpy as np
import uuid

import exceptions
import x3dna
import structure
import basepair
import transform
import util
import io
import motif_type
import settings
import basic_io
import secondary_structure
import secondary_structure_factory as ssf
import rna_structure
import residue
import bead
import motif_state
from primitives import aligner

class Motif(rna_structure.RNAStructure):
    """
    The basic unit of this project stores the 3D coordinates of a RNA Motif
    as well as the 3DNA parameters such as reference frame and origin for
    each basepair

    :param mdir: the path to a motif directory that contains required files
    :type mdir: str

    :param pdb: the path to a pdb file to create this motif object, will also
        creates ref_frames.dat file and dssr out file in current directory
    :type pdb: str

    :param mtype: the enum motif type that this motif is, default UNKNOWN
    :type mtype: motif_type enum

    .. code-block:: python
        #creation from motif dir (recommended)

        #creation from a pdb, generates x3dna files at runtime
        >>> Motif(pdb="test.pdb")
        <Motif(name='test', ends='0')>

    Attributes
    ----------
    `base_pairs` : List of Basepair objects
        All the basepair info determined from 3DNA
    `beads` : List of Bead objects
        All the beads in the 3 bead residue model for all the residues in
        structure object
    `dir` : str
        full path to directory
    `ends` : List of the Basepair objects
        that are at the end of Motif this is critical to assemble motifs
        together its not necessarily the first and last Basepairs
    `structure` : Structure object
        holds 3D coordinate data
    `name` : str
        the name of the directory without entire path
    """

    __slots__ = [
        "_structure",
        "_basepairs",
        "_ends",
        "_name",
        "_score",
        "_protein_beads",
        "_end_ids",
        "_block_end_add",
        "_dot_bracket",
        "_uuid",
        "_mtype",
        "_score"
    ]

    def __init__(self, structure, basepairs, ends, end_ids, name, mtype, score,
                 dot_bracket, block_end_add=0, protein_beads=None, m_uuid=None):

        super(self.__class__, self).__init__(structure, basepairs, ends,
                                             end_ids, name, dot_bracket,
                                             block_end_add, protein_beads)

        self._mtype = mtype
        self._score = score
        self._uuid = m_uuid

        if self._uuid is None:
            self._uuid = uuid.uuid1()

    @classmethod
    def from_str(cls, s, rts):
        spl = s.split("&")
        name = spl[0]
        score = float(spl[1])
        block_end_add = int(spl[2])
        mtype = int(spl[3])
        struc = structure.Structure.from_str(spl[4], rts)
        bp_strs = spl[5].split("@")
        bps = []
        for bp_str in bp_strs[:-1]:
            bps.append(rna_structure.bp_from_str(struc, bp_str))
        ends = []
        end_strs = spl[6].split("@")
        for end_str in end_strs[:-1]:
            ends.append(rna_structure.bp_from_str(struc, end_str))
        end_ids = spl[7].split()
        bead_strs = spl[8].split(";")
        protein_beads = []
        for bead_str in bead_strs[:-1]:
            protein_beads.append(bead.Bead.from_str(bead_str))
        return cls(struc, bps, ends, end_ids, name, mtype, score,
                   block_end_add, spl[7], protein_beads)

    @classmethod
    def copy(cls, m, new_uuid=0):
        s = structure.Structure.copy(m._structure, new_uuid)
        basepairs = []
        ends = []
        protein_beads = [ residue.Bead.copy(b) for b in m._protein_beads ]
        m_uuid = m._uuid

        for bp in m._basepairs:
            if new_uuid:
                bp_res = m.get_bp_res(bp)
                r_pos_1 = m.get_res_index(bp_res[0])
                r_pos_2 = m.get_res_index(bp_res[1])
                res1 = s.get_residue(index=r_pos_1)
                res2 = s.get_residue(index=r_pos_2)
                new_bp = basepair.Basepair.copy_with_new_uuids(bp, res1.uuid, res2.uuid)
                basepairs.append(new_bp)
            else:
                basepairs.append(basepair.Basepair.copy(bp))

        for end in m._ends:
            if new_uuid:
                bp_res = m.get_bp_res(end)
                r_pos_1 = m.get_res_index(bp_res[0])
                r_pos_2 = m.get_res_index(bp_res[1])
                res1 = s.get_residue(index=r_pos_1)
                res2 = s.get_residue(index=r_pos_2)
                new_end = basepair.Basepair.copy_with_new_uuids(end, res1.uuid, res2.uuid)
                ends.append(new_end)
            else:
                ends.append(basepair.Basepair.copy(end))

        if new_uuid:
            m_uuid = uuid.uuid1()

        return cls(s, basepairs, ends, m._end_ids, m._name, m._mtype, m._score,
                   m._block_end_add, m._dot_bracket, protein_beads, m_uuid)

    @classmethod
    def altered_copy(cls, m, name=None, mtype=None):
        if name is None:
            name = m._name
        if mtype is None:
            mtype = m._mtype

        s = structure.Structure.copy(m._structure, new_uuid)
        basepairs = []
        ends = []
        protein_beads = [ residue.Bead.copy(b) for b in m._protein_beads ]

        for bp in m._basepairs:
            if new_uuid:
                bp_res = m.get_bp_res(bp)
                res1 = s.get_residue(num=bp_res[0].num, chain_id=bp_res[0].chain_id,
                                     i_code=bp_res[0].i_code)
                res2 = s.get_residue(num=bp_res[1].num, chain_id=bp_res[1].chain_id,
                                     i_code=bp_res[1].i_code)
                bp = basepair.Basepair.copy_with_new_uuids(bp, res1.uuid, res2.uuid)
                basepairs.append(bp)
            else:
                basepairs.append(basepair.Basepair.copy(bp))

        for end in m._ends:
            i = m._basepairs.index(end)
            ends.append(basepairs[i])

        return cls(s, basepairs, ends, m._end_ids, name, mtype, m._score,
                   m._block_end_add, m._dot_bracket, protein_beads, m._uuid)

    def __repr__(self):
        """
        is called when motif is printed
        """
        return "<Motif(\n\tstructure='%s', \n\tends='%s')>" % (
        self._structure,len(self._ends))

    def to_str(self):
        """
        stringifies motif object
        """
        s  = self._name + "&" + str(self._score) + "&"
        s += str(self._block_end_add) + "&"
        s += str(self._mtype) + "&" + self._structure.to_str() + "&"
        for bp in self._basepairs:
            res1, res2 = self.get_bp_res(bp)
            s += bp.to_str() + ";"
            s += str(res1.num) + "|" + res1.chain_id + "|" + res1.i_code + ";"
            s += str(res2.num) + "|" + res2.chain_id + "|" + res2.i_code + "@"
        s += "&"
        for end in self._ends:
            res1, res2 = self.get_bp_res(end)
            s += end.to_str() + ";"
            s += str(res1.num) + "|" + res1.chain_id + "|" + res1.i_code + ";"
            s += str(res2.num) + "|" + res2.chain_id + "|" + res2.i_code + "@"
        s += "&"
        for end_id in self._end_ids:
            s += end_id + " "
        s += "&"
        s += basic_io.beads_to_str(self._protein_beads)
        s += "&"
        s += self._dot_bracket
        s += "&"
        return s

    def get_state(self):
        basepairs = []
        for bp in self._basepairs:
            sugars = [ np.copy(bp.sugars[0]), np.copy(bp.sugars[1])]
            bp_state = motif_state.Basepair(bp.res1_uuid, bp.res2_uuid, np.copy(bp.r),
                                            np.copy(bp.d), sugars, bp.name, bp.x3dna_bp_type,
                                            bp.bp_type, bp.uuid)
            basepairs.append(bp_state)

        ends = []
        for end in self._ends:
            sugars = [np.copy(end.sugars[0]), np.copy(end.sugars[1])]
            bp_state = motif_state.Basepair(end.res1_uuid, end.res2_uuid, np.copy(end.r),
                                            np.copy(end.d), sugars, end.name, end.x3dna_bp_type,
                                            end.bp_type, end.uuid)
            ends.append(bp_state)

        # keep only watson and crick bps
        #bps = []
        #for bp in basepairs:
        #    if bp.bp_type == "cW-W":
        #        bps.append(bp)
        ms = motif_state.Motif(self._structure.get_state(), basepairs, ends, self._end_ids,
                               self._name, self._mtype, self._score, self._dot_bracket,
                               self._block_end_add, self._uuid)
        return ms

    @property
    def score(self):
        return self._score

    @property
    def mtype(self):
        return self._mtype

    @property
    def uuid(self):
        return self._uuid


class MotifAligner(aligner.Aligner):
    __slots__ = [
        '_r',
        '_d',
        '_t',
        '_m_end',
        '_bp_pos_diff',
        '_dist_1',
        '_dist_2',
        '_sugar_diff_1',
        '_sugar_diff_2'
    ]

    def __init__(self):
        super(self.__class__, self).__init__()
        self._bp_pos_diff = 0
        self._m_end = None
        self._dist_1 = 0
        self._dist_2 = 0
        self._sugar_diff_1 = 0
        self._sugar_diff_2 = 0

    def get_aligned_motif(self, ref_bp, m):
        m_copy = Motif.copy(m)
        self.align(ref_bp, m_copy)
        return m_copy

    def align(self, ref_bp, m):
        self._m_end = m.get_end(0)

        # calculate rotation before ref_bp (where we are going) and current
        # base pair (where we are)
        self._r = util.unitarize(ref_bp.r.T.dot(self._m_end.r))
        self._d = -self._m_end.d
        self._t = transform.Transform(self._r, self._d)
        m.transform(self._t)
        self._bp_pos_diff = ref_bp.d - self._m_end.d
        m.move(self._bp_pos_diff)

        # alignment is by center of basepair, it can be slightly improved by
        # aligning the c1' sugars
        self._dist_1 = util.distance(self._m_end.res1_sugar, ref_bp.res1_sugar)
        self._dist_2 = util.distance(self._m_end.res2_sugar, ref_bp.res1_sugar)

        if self._dist_1 < self._dist_2:
            self._sugar_diff_1 = ref_bp.res1_sugar - self._m_end.res1_sugar
            self._sugar_diff_2 = ref_bp.res2_sugar - self._m_end.res2_sugar
        else:
            self._sugar_diff_1 = ref_bp.res1_sugar - self._m_end.res2_sugar
            self._sugar_diff_2 = ref_bp.res2_sugar - self._m_end.res1_sugar

        if self._dist_1 < 5 or self._dist_2 < 5:
            m.move((self._sugar_diff_1 + self._sugar_diff_2) / 2)


def align_motif(ref_bp_state, motif_end, motif):
    """
    This is the workhorse of the entire suite. Aligns one end of a motif to
    the reference frame and origin of a Basepair object.

    :param ref_bp: the base pair that the motif end is going to align too
    :param motif_end: the motif end basepair to overly with the ref_bp
    :param motif: the motif object that you want to align

    :type ref_bp: Basepair object
    :type motif_end: Basepair object
    :type motif: Motif object
    """

    r1 , r2 = ref_bp_state.r , motif_end.r
    r = util.unitarize(r1.T.dot(r2))
    trans = -motif_end.d
    t = transform.Transform(r, trans)
    motif.transform(t)
    bp_pos_diff = ref_bp_state.d - motif_end.d
    motif.move(bp_pos_diff)

    #alignment is by center of basepair, it can be slightly improved by
    #aligning the c1' sugars
    res1_coord, res2_coord = motif_end.sugars
    ref_res1_coord, ref_res2_coord = ref_bp_state.sugars

    dist1 = util.distance(res1_coord, ref_res1_coord)
    dist2 = util.distance(res2_coord, ref_res1_coord)

    if dist1 < dist2:
        sugar_diff_1 = ref_res1_coord - res1_coord
        sugar_diff_2 = ref_res2_coord - res2_coord
    else:
        sugar_diff_1 = ref_res1_coord - res2_coord
        sugar_diff_2 = ref_res2_coord - res1_coord

    if dist1 < 5 or dist2 < 5:
        motif.move( (sugar_diff_1 + sugar_diff_2) / 2 )


def get_aligned_motif(ref_bp, motif_end, m, sterics=1):

    motif_end_index = m.get_end_index(motif_end.name)
    m_copy = Motif.copy(m)
    motif_end = m_copy.get_end(motif_end_index)
    align_motif(ref_bp, motif_end, m_copy)

    return m_copy


def clash_between_motifs(m1, m2, clash_radius=settings.CLASH_RADIUS):
    for r1 in m1:
        for b1 in r1.iter_beads():
            for r2 in m2:
                for b2 in r2.iter_beads():
                    if b1.btype == residue.BeadType.PHOS or \
                       b2.btype == residue.BeadType.PHOS:
                        continue
                    dist = util.distance(b1.center, b2.center)
                    if dist < clash_radius:
                        return 1
    return 0

