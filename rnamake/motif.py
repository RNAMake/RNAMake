import itertools
import numpy as np
import uuid

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

    def __init__(self, r_struct=None):
        self.beads, self.score, self.mtype, self.basepairs = [], 0, motif_type.UNKNOWN, []
        self.path, self.name, self.ends = "", "", []
        self.end_ids = []
        self.structure = structure.Structure()
        self.secondary_structure = secondary_structure.Motif()
        self.block_end_add = 0
        self.id = uuid.uuid1()
        self.protein_beads = []

        if r_struct is not None:
            self.__dict__.update(r_struct.__dict__)

    def __repr__(self):
        """
        is called when motif is printed
        """
        return "<Motif(\n\tstructure='%s', \n\tends='%s')>" % (
        self.structure,len(self.ends))

    def to_str(self):
        """
        stringifies motif object
        """
        s = self.path + "&" + self.name + "&" + str(self.score) + "&" +  \
            str(self.block_end_add) + "&" + \
            str(self.mtype) + "&" + self.structure.to_str() + "&"
        for bp in self.basepairs:
            s += bp.to_str() + "@"
        s += "&"
        for end in self.ends:
            index = self.basepairs.index(end)
            s += str(index) + " "
        s += "&"
        for end_id in self.end_ids:
            s += end_id + " "
        s += "&"
        s += self.secondary_structure.to_str()
        s += "&"
        return s

    def transform(self, t):
        """
        perform an transformation of both structure and basepairs
        """
        r_T = t.rotation().T
        for bp in self.basepairs:
            transformed = np.dot(bp.state().r, r_T)
            bp.state().r = transformed

        self.structure.transform(t)

    def move(self, p):
        """
        wrapper for self.structure.move
        """
        return self.structure.move(p)

    def copy(self):
        """
        performs a deep copy of this motif
        """
        cmotif = Motif()
        cmotif.name      = self.name
        cmotif.path      = self.path
        cmotif.score     = self.score
        cmotif.mtype     = self.mtype
        cmotif.id        = self.id
        cmotif.structure = self.structure.copy()
        cmotif.beads     = [b.copy() for b in self.beads]
        cmotif.end_ids   = list(self.end_ids)
        cmotif.secondary_structure = self.secondary_structure.copy()
        cmotif.block_end_add = self.block_end_add

        for bp in self.basepairs:
            new_res1 = cmotif.get_residue(uuid=bp.res1.uuid)
            new_res2 = cmotif.get_residue(uuid=bp.res2.uuid)
            # hopefully this doesnt happen anymore
            if new_res1 is None or new_res2 is None:
                raise ValueError("could not find a residue during copy")
            new_r = np.copy(bp.bp_state.r)
            new_bp = basepair.Basepair(new_res1, new_res2, new_r, bp.bp_type)
            new_bp.uuid = bp.uuid
            cmotif.basepairs.append(new_bp)

        for end in self.ends:
            index = self.basepairs.index(end)
            cmotif.ends.append(cmotif.basepairs[index])

        return cmotif

    def get_state(self):
        end_res = self.ends[0].residues()
        res_states = []
        for r in self.residues():
            if r in end_res:
                res_states.append(residue.get_residue_state(r,
                                                    residue.ResidueStateType.END))
            else:
                res_states.append(residue.get_residue_state(r))
        beads = []
        for r in res_states:
            if r.type != residue.ResidueStateType.END:
                beads.extend(r.beads)

        ends = [None for x in self.ends]
        end_names = []
        for i, end in enumerate(self.ends):
            ends[i] = end.state()
            end_names.append(end.name())
        return MotifState(self.name, end_names, self.end_ids, ends, beads,
                          self.score, len(self.residues()), self.block_end_add,
                          res_states)

    def end_index_with_id(self, id):
        for i, end_id in enumerate(self.end_ids):
            if i == self.block_end_add:
                continue
            if id == end_id:
                return i

        raise ValueError("no end id: " + id + " in motif: " + self.name)

    # TODO rename since it updates all unique indenfiers
    def new_res_uuids(self):
        self.id = uuid.uuid1()
        for i, r in enumerate(self.residues()):
            ss_r = self.secondary_structure.get_residue(uuid=r.uuid)
            r.new_uuid()
            ss_r.uuid = r.uuid
        for bp in self.basepairs:
            bp.uuid = uuid.uuid1()

    def to_secondary_structure(self):
        ss = ssf.factory.secondary_structure_from_motif(self)
        return ss

class MotifState(object):
    __slots__ = ['name', 'end_names', 'end_ids', 'end_states',
                 'beads', 'score', 'size', 'block_end_add', 'residues']

    def __init__(self, name, end_names, end_ids, end_states,
                 beads, score, size, block_end_add, residues):
        self.name, self.end_states, self.beads = name, end_states, beads
        self.score, self.size = score, size
        self.end_names, self.end_ids = end_names, end_ids
        self.block_end_add = block_end_add
        self.residues = residues

    def to_str(self):
        s = self.name + "|" + str(self.score) + "|" + str(self.size) + "|"
        s += str(self.block_end_add) + "|"
        s += basic_io.points_to_str(self.beads) + "|"
        #for r in self.residues:
        #    s += r.to_str()
        #    s += "E"
        #s += "|"
        s += ",".join(self.end_names) + "|"
        s += ",".join(self.end_ids) + "|"
        for state in self.end_states:
            s += state.to_str() + "|"
        return s

    def get_end_state(self, name):
        for i, n in enumerate(self.end_names):
            if i == self.block_end_add:
                continue
            if n == name:
                return self.end_states[i]
        raise ValueError("cannot get end state of name: " + name)

    def end_index_with_id(self, id):
        for i, end_id in enumerate(self.end_ids):
            if i == self.block_end_add:
                continue
            if id == end_id:
                return i

        raise ValueError("no matching end id in motif")

    def copy(self):
        end_states = [end.copy() for end in self.end_states]
        residues = [r.copy() for r in self.residues]
        beads = []
        for r in residues:
            if r.type != residue.ResidueStateType.END:
                beads.extend(r.beads)
        #print "beads", len(beads)
        return MotifState(self.name, self.end_names, self.end_ids, end_states,
                          self.beads, self.score, self.size, self.block_end_add, residues)

    def update_res_uuids(self, res):
        for i, r in enumerate(res):
            self.residues[i].uuid = r.uuid


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
    return m


def str_to_motif_state(s):
    spl = s.split("|")
    name, score, size = spl[0], float(spl[1]), float(spl[2])
    block_end_add = int(spl[3])
    residues = []
    #beads = basic_io.str_to_points(spl[4])
    res_str = spl[4].split("E")
    beads = []
    for r_str in res_str[:-1]:
        r_spl = r_str.split(",")
        points = basic_io.str_to_points(r_spl[0])
        r = residue.ResidueState2Bead(None, points[0], points[1], int(r_spl[1]))
        residues.append(r)
        if r.type != residue.ResidueStateType.END:
            beads.extend(r.beads)
    end_names = spl[5].split(",")
    end_ids = spl[6].split(",")
    end_states = []
    for i in range(7, len(spl)-1):
            end_states.append(basepair.str_to_basepairstate(spl[i]))

    return MotifState(name, end_names, end_ids, end_states, beads, score, size, block_end_add, residues)


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

    r1 , r2 = ref_bp_state.r , motif_end.state().r
    r = util.unitarize(r1.T.dot(r2))
    trans = -motif_end.state().d
    t = transform.Transform(r, trans)
    motif.transform(t)
    bp_pos_diff = ref_bp_state.d - motif_end.state().d
    motif.move(bp_pos_diff)

    #alignment is by center of basepair, it can be slightly improved by
    #aligning the c1' sugars
    res1_coord, res2_coord = motif_end.c1_prime_coords()
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

    if sterics:
        motif.get_beads([motif_end])


def get_aligned_motif(ref_bp, motif_end, motif, sterics=1):

    motif_end_index = motif.ends.index(motif_end)
    m_copy = motif.copy()
    motif_end = m_copy.ends[motif_end_index]

    align_motif(ref_bp.state(), motif_end, m_copy)

    return m_copy


def align_motif_state(ref_bp_state, org_state):
    r, t = ref_bp_state.get_transforming_r_and_t_w_state(org_state.end_states[0])
    t += ref_bp_state.d

    for i, s in enumerate(org_state.end_states):
        new_r, new_d, new_sug = s.get_transformed_state(r, t)
        org_state.end_states[i].set(new_r,new_d,new_sug)


def get_aligned_motif_state(ref_bp_state, cur_state, org_state):
    r, t = ref_bp_state.get_transforming_r_and_t_w_state(org_state.end_states[0])
    t += ref_bp_state.d

    for i, s in enumerate(org_state.end_states):
        new_r, new_d, new_sug = s.get_transformed_state(r, t)
        cur_state.end_states[i].set(new_r, new_d, new_sug)

    """cur_state.beads = []
    r_T = r.T
    for r in cur_state.residues:
        r.update(r_T, t)
        if r.type != residue.ResidueStateType.END:
            cur_state.beads.extend(r.beads)
    """
    if len(org_state.beads) > 0:
        cur_state.beads = np.dot(org_state.beads, r.T) + t



