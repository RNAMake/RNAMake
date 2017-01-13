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
                 block_end_add=0, dot_bracket=None, protein_beads=None, m_uuid=None):

        super(self.__class__, self).__init__(structure, basepairs, ends,
                                             end_ids, name, block_end_add,
                                             dot_bracket, protein_beads)

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
        ends = [ bps[int(i)] for i in spl[6].split() ]
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
                bp = basepair.Basepair.copy_with_new_uuids(bp, res1.uuid, res2.uuid)
                basepairs.append(bp)
            else:
                basepairs.append(basepair.Basepair.copy(bp))

        for end in m._ends:
            i = m._basepairs.index(end)
            ends.append(basepairs[i])

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
            index = self._basepairs.index(end)
            s += str(index) + " "
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
                                            np.copy(bp.d), sugars, bp.name, bp.bp_type,
                                            bp.uuid)
            basepairs.append(bp_state)

        ends = []
        for end in self._ends:
            i = self._basepairs.index(end)
            ends.append(basepairs[i])

        # keep only watson and crick bps
        bps = []
        for bp in basepairs:
            if bp.bp_type == "cW-W":
                bps.append(bp)

        ms = motif_state.Motif(self._structure.get_state(), bps, ends, self._end_ids,
                               self._name, self._mtype, self._score,
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


def file_to_motif(path):
    try:
        f = open(path)
        l = f.readline()
        f.close()
    except IOError:
        raise IOError("cannot find path to open motif from")

    return str_to_motif(l)


def str_to_motif(s):
    """
    creates motif from stringified motif, this is created by motif.to_str()

    :param s: stringified motif
    :type s: str
    """
    spl = s.split("&")
    m = Motif()
    m.path = spl[0]
    m.name = spl[1]
    m.score = float(spl[2])
    m.block_end_add = int(spl[3])
    m.mtype = int(spl[4])
    m.structure = io.str_to_structure(spl[5])
    m.basepairs = []
    m.id = uuid.uuid1()

    basepair_str = spl[6].split("@")
    for bp_str in basepair_str[:-1]:
        bp_spl = bp_str.split(",")
        res_spl = bp_spl[0].split("-")
        res1_id, res1_num = res_spl[0][0], int(res_spl[0][1:])
        res2_id, res2_num = res_spl[1][0], int(res_spl[1][1:])
        res1 = m.get_residue(num=res1_num, chain_id=res1_id)
        res2 = m.get_residue(num=res2_num, chain_id=res2_id)
        state = basepair.str_to_basepairstate(bp_spl[1])
        bp = basepair.Basepair(res1, res2, state.r, bp_spl[2])
        m.basepairs.append(bp)

    end_indexes = spl[7].split()
    for index in end_indexes:
        m.ends.append(m.basepairs[int(index)])
    end_ids = spl[8].split()
    m.end_ids = end_ids
    m.secondary_structure = secondary_structure.str_to_motif(spl[9])
    ss_res = m.secondary_structure.residues()
    for i, r in enumerate(m.residues()):
        ss_res[i].uuid = r.uuid
    for b_str in spl[10].split(";"):
        b_spl = b_str.split(",")
        if len(b_spl) < 2:
            continue
        b = residue.Bead(basic_io.str_to_point(b_spl[0]), int(b_spl[1]))
        m.protein_beads.append(b)
    return m


def str_to_motif_state(s):
    spl = s.split("|")
    name, score, size = spl[0], float(spl[1]), float(spl[2])
    block_end_add = int(spl[3])
    residues = []
    beads = basic_io.str_to_points(spl[4])
    end_names = spl[5].split(",")
    end_ids = spl[6].split(",")
    end_states = []
    for i in range(7, len(spl)-1):
            end_states.append(basepair.str_to_basepairstate(spl[i]))

    return MotifState(name, end_names, end_ids, end_states, beads, score, size, block_end_add)


def align_motif(ref_bp_state, motif_end, motif, sterics=1):
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
    for r1 in m1.iter_res():
        for b1 in r1.iter_beads():
            for r2 in m2.iter_res():
                for b2 in r2.iter_beads():
                    if b1.btype == residue.BeadType.PHOS or \
                       b2.btype == residue.BeadType.PHOS:
                        continue
                    dist = util.distance(b1.center, b2.center)
                    if dist < clash_radius:
                        return 1
    return 0

