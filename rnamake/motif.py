import itertools
import x3dna
import structure
import basepair
import transform
import util
import io
import motif_type
import settings
import basic_io
import numpy as np

class Motif(object):
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
        >>> Motif(pdb=test.pdb")
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

    def __init__(self):
        self.beads, self.score, self.mtype, self.basepairs = [], 0, motif_type.UNKNOWN, []
        self.path, self.name, self.ends = "", "", []
        self.ss_chains = []
        self.end_ids = []
        self.structure = structure.Structure()

    def __repr__(self):
        """
        is called when motif is printed
        """
        return "<Motif(\n\tstructure='%s', \n\tends='%s')>" % (
        self.structure,len(self.ends))

    def get_basepair(self, bp_uuid=None, res1=None, res2=None, uuid1=None,
                     uuid2=None, name=None):
        """
        locates a Basepair object based on residue objects or uuids if nothing
        is supplied you will get back all the basepairs in the motif. The way
        to make sure you will only get one basepair back is to supply BOTH
        res1 and res2 OR uuid1 and uuid2, I have left it open like this
        because it is sometimes useful to get all basepairs that a single
        residue is involved

        :param res1: First residue
        :param res2: Second residue
        :param uuid1: First residue uuid
        :param uuid2: Second residue uuid

        :type res1: Residue object
        :type res2: Residue object
        :type uuid1: uuid object
        :type uuid2: uuid object
        """
        alt_name = None
        if name:
            name_spl = name.split("-")
            alt_name = name_spl[1] + "-" + name_spl[0]

        found = []
        for bp in self.basepairs:
            if bp_uuid is not None and bp_uuid != bp.uuid:
                continue
            if res1 is not None and (res1 != bp.res1 and res1 != bp.res2):
                continue
            if res2 is not None and (res2 != bp.res1 and res2 != bp.res2):
                continue
            if uuid1 is not None and \
               (uuid1 != bp.res1.uuid and uuid1 != bp.res2.uuid):
                continue
            if uuid2 is not None and \
               (uuid2 != bp.res1.uuid and uuid2 != bp.res2.uuid):
                continue
            if name is not None and \
               (name != bp.name() and alt_name != bp.name()):
                continue
            found.append(bp)
        return found

    def get_beads(self, excluded_ends=None, excluded_res=None):
        excluded = []
        if excluded_ends:
            for end in excluded_ends:
                excluded.extend(end.residues())

        if excluded_res:
            excluded.extend(excluded_res)

        self.beads = self.structure.get_beads(excluded)
        return self.beads

    def sequence(self):
        seqs = [x.seq for x in self.ss_chains]
        return "&".join(seqs)

    def secondary_structure(self):
        sss = [x.ss for x in self.ss_chains]
        return "&".join(sss)

    def to_str(self):
        """
        stringifies motif object
        """
        s = self.path + "&" + self.name + "&" + str(self.score) + "&" + \
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
        return s

    def to_pdb_str(self):
        """
        returns pdb formatted string of motif's structure object
        """
        return self.structure.to_pdb_str()

    def to_pdb(self, fname="motif.pdb"):
        """
        writes the current motif's structure to a pdb
        """
        return self.structure.to_pdb(fname)

    def get_residue(self, num=None, chain_id=None, i_code=None, uuid=None):
        """
        wrapper to self.structure.get_residue()
        """
        return self.structure.get_residue(num=num, chain_id=chain_id,
                                          i_code=i_code, uuid=uuid)

    def residues(self):
        """
        wrapper to self.structure.residues()
        """
        return self.structure.residues()

    def chains(self):
        return self.structure.chains

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
        cmotif.structure = self.structure.copy()
        cmotif.beads     = [b.copy() for b in self.beads]
        cmotif.ss_chains = list(self.ss_chains)
        cmotif.end_ids   = list(self.end_ids)

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
        beads = self.get_beads([self.ends[0]])
        bead_centers = []
        for b in beads:
            if b.btype == 0:
                continue
            bead_centers.append(b.center)
        ends = [None for x in self.ends]
        for i, end in enumerate(self.ends):
            ends[i] = end.state()
        return MotifState(self.name, ends, bead_centers, self.score, len(self.residues()))


class MotifState(object):
    __slots__ = ['name', 'end_states', 'beads', 'score', 'size']

    def __init__(self, name, end_states, beads, score, size):
        self.name, self.end_states, self.beads = name, end_states, beads
        self.score, self.size = score, size

    def to_str(self):
        s = self.name + "|" + str(self.score) + "|" + str(self.size) + "|"
        s += basic_io.points_to_str(self.beads) + "|"
        for state in self.end_states:
            s += state.to_str() + "|"
        return s

    def copy(self):
        end_states = [end.copy() for end in self.end_states]
        beads = np.copy(self.beads)

        return MotifState(self.name, end_states, beads, self.score, self.size)


class MotifArray(object):
    __slots__ = ['motifs']

    def __init__(self, motifs=[]):
        self.motifs = motifs

    def add(self, m):
        self.motifs.append(m)

    def to_str(self):
        s = ""
        for m in self.motifs:
            s += m.to_str() + "$"
        return s


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
    m.mtype = int(spl[3])
    m.structure = io.str_to_structure(spl[4])
    m.basepairs = []

    basepair_str = spl[5].split("@")
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

    end_indexes = spl[6].split()
    for index in end_indexes:
        m.ends.append(m.basepairs[int(index)])
    end_ids = spl[7].split()
    m.end_ids = end_ids

    return m


def str_to_motif_state(s):
    spl = s.split("|")
    name, score, size = spl[0], float(spl[1]), float(spl[2])
    beads = basic_io.str_to_points(spl[3])
    end_states = []
    for i in range(4, len(spl)-1):
            end_states.append(basepair.str_to_basepairstate(spl[i]))

    return MotifState(name, end_states, beads, score, size)


def str_to_motif_array(str):
    spl = str.split("$")
    motifs = []
    for s in spl:
        motifs.append(str_to_motif(s))
    return MotifArray(motifs)


def align_motif(ref_bp, motif_end, motif, sterics=1):
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

    r1 , r2 = ref_bp.state().r , motif_end.state().r
    r = util.unitarize(r1.T.dot(r2))
    trans = -motif_end.state().d
    t = transform.Transform(r, trans)
    motif.transform(t)
    bp_pos_diff = ref_bp.state().d - motif_end.state().d
    motif.move(bp_pos_diff)

    #alignment is by center of basepair, it can be slightly improved by
    #aligning the c1' sugars
    res1_coord, res2_coord = motif_end.c1_prime_coords()
    ref_res1_coord, ref_res2_coord = ref_bp.c1_prime_coords()

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

    align_motif(ref_bp, motif_end, m_copy)

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
        cur_state.end_states[i].set(new_r,new_d,new_sug)

    cur_state.beads = np.dot(cur_state.beads, r.T) + t















