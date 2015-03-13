import itertools
import x3dna
import structure
import basepair
import transform
import util
import io
import motif_type
import settings
import numpy as np

class MotifException(Exception):
    pass


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

    def __init__(self, mdir=None, pdb=None, mtype=motif_type.UNKNOWN):
        if mdir is not None and pdb is not None:
            raise ValueError("cannot initiate a Motif with both a mdir and" +
                             "a pdb")

        self.beads, self.score, self.mtype, self.basepairs = [], 0, mtype, []
        self.mdir, self.name, self.ends = "", "", []
        self.cached_rotations = []
        self._setup(mdir, pdb)

    def __repr__(self):
        """
        is called when motif is printed
        """
        return "<Motif(name='%s', ends='%s')>" % (
        self.name,len(self.ends))

    def _setup(self, mdir=None, pdb=None):
        #nothing to do
        if mdir is None and pdb is None:
            self.structure = structure.Structure()
            return
        # supplied a motif directory that already contains ref_frames.dat and
        # dssr output file
        if mdir:
            filename = util.filename(mdir)
            self.mdir = mdir
            self.name = filename
            self.structure = structure.Structure(mdir + "/" + filename + ".pdb")
        # supplied only a pdb, have to create the ref_frames and dssr output
        # locally
        if pdb:
            mdir = util.base_dir(pdb)
            filename = util.filename(pdb)
            self.mdir = mdir
            self.name = filename[:-4]
            self.structure = structure.Structure(pdb)

        self.basepairs = self._setup_basepairs()
        self.setup_basepair_ends()
        self._cache_basepair_frames()

    def _setup_basepairs(self):
        """
        gets x3dna data on basepairing information and then interwines it
        with the structural information stored in structure for simpler
        retreival of data
        """
        x3dna_parser = x3dna.X3dna()
        mdir = self.mdir + "/"
        if mdir == "//":
            mdir = ""
        x_basepairs = x3dna_parser.get_basepairs(mdir + self.name)
        basepairs = []
        for xbp in x_basepairs:
            res1 = self.structure.get_residue(num=xbp.res1.num,
                                              chain_id=xbp.res1.chain_id,
                                              i_code=xbp.res1.i_code)

            res2 = self.structure.get_residue(num=xbp.res2.num,
                                              chain_id=xbp.res2.chain_id,
                                              i_code=xbp.res2.i_code)

            if res1 is None or res2 is None:
                raise ValueError("cannot find residues in basepair")

            bp = basepair.Basepair(res1, res2, xbp.r, xbp.bp_type)
            self._assign_bp_primes(bp)
            basepairs.append(bp)

        return basepairs

    def _assign_bp_primes(self, bp):
        """
        This is legacy code not sure if I still need it The purpose is to
        determine which residue in a Basepair object is going in the 5' and 3'
        direction. How this is calculated is position from a chain end. If
        residue 1 is positioned in the half the chain its called the 5' end and
        if its in the second half its 3'. If its exactly at the 1/2 point then
        the second residue is checked in the same manner.

        :param bp: an end basepair
        :type bp: Basepair object
        """
        res1_pos, res2_pos, res1_total, res2_total = 0, 0, None, None
        for c in self.structure.chains:
            for i, r in enumerate(c.residues):
                if bp.res1 == r:
                    res1_pos = i
                    res1_total = len(c.residues)
                elif bp.res2 == r:
                    res2_pos = i
                    res2_total = len(c.residues)

        if res1_pos > res1_total/2 or res2_pos < res2_total/2:
            bp.res1, bp.res2 = bp.res2, bp.res1

    def _cache_basepair_frames(self):
        self.cached_rotations = []
        for bp in self.basepairs:
            self.cached_rotations.append(np.copy(bp.state().r))

    def setup_basepair_ends(self):
        # TODO revisit this code, is it really necessary?
        # get chain ends

        chain_ends = []
        for c in self.structure.chains:
            chain_ends.append(c.first())
            if len(c) > 1:
                chain_ends.append(c.last())

        seen_res_bps = {res: [] for res in chain_ends}
        for bp in self.basepairs:
            for ce1 in chain_ends:
                for ce2 in chain_ends:
                    if bp.res1 == ce1 and bp.res2 == ce2:
                        seen_res_bps[ce1].append(bp)
                        seen_res_bps[ce2].append(bp)

        lists = filter(lambda l: len(l) > 0, seen_res_bps.values())
        combos = itertools.product(*lists)
        best, best_count = [], 0
        for c in combos:
            bps, res = [], []
            for bp in c:
                if bp in bps:
                    continue
                fail = 0
                for r in bp.residues():
                    if r in res:
                        fail = 1
                        break
                    res.append(r)
                if fail:
                    continue
                bps.append(bp)
                if len(res) > best_count:
                    best_count = len(res)
                    best = bps

        # for reproducability
        best.sort(key=lambda x: x.name())
        self.ends = best
        return best

    def get_basepair(self, bp_uuid=None, res1=None, res2=None, uuid1=None,
                     uuid2=None):
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
        seq = ""
        for c in self.chains():
            for r in c.residues:
                seq += r.rtype.name[0]
            seq += "&"
        return seq[:-1]

    def secondary_structure(self):
        structure = ""
        seen_res = {}
        seen_bp = {}
        saved_bp = None
        for c in self.chains():
            for r in c.residues:
                ss = ""
                bps = self.get_basepair(res1=r)
                is_bp = 0
                for bp in bps:
                    partner_res = bp.partner(r)
                    is_bp = 1
                    passes = 0
                    saved_bp = None
                    if util.wc_bp(bp) and bp.bp_type == "cW-W":
                        passes = 1
                    if util.gu_bp(bp) and bp.bp_type == "cW-W":
                        passes = 1

                    if passes:
                        saved_bp = bp
                        if   bp not in seen_bp and r not in seen_res and \
                             partner_res not in seen_res:
                            seen_res[r] = 1
                            ss = "("
                        elif partner_res in seen_res:
                            if seen_res[partner_res] > 1:
                                ss = "."
                            else:
                                ss = ")"
                                seen_res[r] = 1
                                seen_res[partner_res] += 1
                                break
                    elif r not in seen_res:
                        ss = "."

                if not is_bp:
                    ss = "."
                if saved_bp is not None:
                    seen_bp[saved_bp] = 1
                structure += ss
            structure += "&"
        return structure[:-1]

    def to_str(self):
        """
        stringifies motif object
        """
        s = self.mdir + "&" + self.name + "&" + str(self.score) + "&" + \
            str(self.mtype) + "&" + self.structure.to_str() + "&"
        for bp in self.basepairs:
            s += bp.to_str() + "@"
        s += "&"
        for end in self.ends:
            index = self.basepairs.index(end)
            s += str(index) + " "
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
        cmotif.name = self.name
        cmotif.mdir = self.mdir
        cmotif.score = self.score
        cmotif.mtype = self.mtype
        cmotif.structure = self.structure.copy()
        cmotif.beads = [b.copy() for b in self.beads]
        # cmotif.cached_rotations = self.cached_rotations
        for bp in self.basepairs:
            new_res1 = cmotif.get_residue(uuid=bp.res1.uuid)
            new_res2 = cmotif.get_residue(uuid=bp.res2.uuid)
            # hopefully this doesnt happen anymore
            if new_res1 is None or new_res2 is None:
                raise MotifException("could not find a residue during copy")
            new_r = np.copy(bp.bp_state.r)
            new_bp = basepair.Basepair(new_res1, new_res2, new_r, bp.bp_type)
            new_bp.designable = bp.designable
            new_bp.flipped = bp.flipped
            new_bp.uuid = bp.uuid
            cmotif.basepairs.append(new_bp)

        for end in self.ends:
            index = self.basepairs.index(end)
            cmotif.ends.append(cmotif.basepairs[index])
        cmotif._cache_basepair_frames()

        return cmotif

    def reset(self):
        """
        reset both the structure and basepair rotations, so a new
        transformation object can be applied
        """

        for i,bp in enumerate(self.basepairs):
            bp.state().r = self.cached_rotations[i]

        for end in self.ends:
            end.flip(0)

        self.structure.restore_coords()
        self.beads = []

        for i, e in enumerate(self.ends):
            if e.uuid == end.uuid:
                return i
        raise ValueError("end is not a end of this motif")


def str_to_motif(s):
    """
    creates motif from stringified motif, this is created by motif.to_str()

    :param s: stringified motif
    :type s: str
    """
    spl = s.split("&")
    m = Motif()
    m.mdir = spl[0]
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
        bp.designable = int(bp_spl[3])
        bp.flipped = int(bp_spl[4])
        m.basepairs.append(bp)

    end_indexes = spl[6].split()
    for index in end_indexes:
        m.ends.append(m.basepairs[int(index)])
    m._cache_basepair_frames()
    return m


def align_motif(ref_bp, motif_end, motif):
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

def ref_motif():
    path = settings.RESOURCES_PATH + "/start"
    m = Motif(path)
    return m
