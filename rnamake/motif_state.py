import uuid
import numpy as np

import primitives.residue
import primitives.chain
import primitives.structure
import primitives.basepair
from primitives.ensemble import Ensemble, EnsembleMember
import primitives.rna_structure
import bead
import exceptions
import basic_io
import util
import rna_structure
import transform

class Residue(primitives.residue.Residue):
    __slots__ = [
        "_name",
        "_num",
        "_chain_id",
        "_i_code",
        "_uuid",
        "_beads"]

    def __init__(self, name, num, chain_id, i_code, beads, r_uuid=None):
        super(self.__class__, self).__init__(name, num, chain_id, i_code, r_uuid)
        self._beads = beads

    def iter_beads(self):
        return self._beads.__iter__()

    @classmethod
    def from_str(cls, s):
        spl = s.split(",")
        beads = []
        for i in range(4, len(spl)-1, 2):
            b = bead.Bead.from_str(spl[i] + "," + spl[i+1])
            beads.append(b)

        return cls(spl[0], int(spl[1]), spl[2], spl[3], beads)

    @classmethod
    def copy(cls, r, new_uuid=0):
        r_uuid = r._uuid
        if new_uuid:
            r_uuid = uuid.uuid1()

        beads = [ bead.Bead.copy(b) for b in r._beads]
        return cls(r._name, r._num, r._chain_id, r._i_code, beads, r_uuid)

    def to_str(self):
        bead_str = ""
        for b in self._beads:
            bead_str += b.to_str() + ","

        return self._name + "," + str(self._num) + "," + \
               str(self._chain_id) + "," + str(self._i_code) + "," + bead_str

    def move(self, p):
        for b in self._beads:
            b.move(p)

    def transform(self, t):
        for b in self._beads:
            b.transform(t)

    def fast_transform(self, r, t):
        for b in self._beads:
            b.fast_transform(r, t)

    def num_beads(self):
        return len(self._beads)


class Chain(primitives.chain.Chain):
    __slots__ = ["_residues"]

    def __init__(self, residues):
        super(self.__class__, self).__init__(residues)

    @classmethod
    def from_str(cls, s):
        spl = s.split(";")
        residues = []
        for r_str in spl[:-1]:
            r = Residue.from_str(r_str)
            residues.append(r)
        return cls(residues)

    @classmethod
    def copy(cls, c, new_uuid=0):

        residues = [Residue.copy(r, new_uuid) for r in c]
        return cls(residues)

    def to_str(self):
        s = ""
        for r in self._residues:
            s += r.to_str() + ";"
        return s

    def move(self, p):
        for r in self._residues:
            r.move(p)

    def transform(self, t):
        for r in self._residues:
            r.transform(t)

    def fast_transform(self, r, t):
        for res in self._residues:
            res.fast_transform(r, t)


class Structure(primitives.structure.Structure):
    __slots__ = [
        "_chains",
        "_residues"
    ]

    def __init__(self, chains):
        super(self.__class__, self).__init__(chains)

    @classmethod
    def from_str(cls, s):
        spl = s.split(":")
        chains = []
        for c_str in spl[:-1]:
            c = Chain.from_str(c_str)
            chains.append(c)

        return cls(chains)

    @classmethod
    def copy(cls, s, new_uuid=0):
        """
        creates a deep copy of this structure

        :returns: copy of Structure object
        :rtype: Structure
        """
        chains = []
        for c in s._chains:
            cc = Chain.copy(c, new_uuid)
            chains.append(cc)

        return cls(chains)

    def to_str(self):
        """
        Stringifes Structure object

        :return: str
        """

        s = ""
        for c in self._chains:
            s += c.to_str() + ":"
        return s

    def move(self, p):
        for c in self._chains:
            c.move(p)

    def transform(self, t):
        for c in self._chains:
            c.transform(t)

    def fast_transform(self, r, t):
        for c in self._chains:
            c.fast_transform(r, t)


class Basepair(primitives.basepair.Basepair):
    """
    :param res1: First residue in basepair
    :param res2: Second residue in basepair
    :param r: Reference frame of basepair
    :param bp_type: X3dna basepair type, default "c..."

    :type res1: residue.Residue
    :type res2: residue.Residue
    :type r: np.array
    :type bp_type: str

    :attributes:

    `res1` : residue.Residue
        First residue in basepair
    `res2` : residue.Residue
        Second residue in basepair
    `bp_type`: str
        X3dna basepair type
    `atoms`: list of atom.Atoms
        All atoms from both res1 and res2. Grouped together for easy manipulation
    `bp_state`: BasepairState
        Contains the orientation, origin and sugar. This is all information
        required to align to this basepair and can be used independently from
        the reset of the basepair.
    `uuid`: uuid.uuid1()
        unique id to indentify this basepair when locating it in a motif or
        pose

    :examples:

    ..  code-block:: python

        # build basepair from stratch
        >>> from rnamake.unittests import instances
        >>> from rnamake import basepair
        >>> import numpy as np
        >>> s = instances.structure()
        >>> b = basepair.Basepair(s.get_residue(num=103),
                                  s.get_residue(num=104),
                                  np.eye(3))
        >>> print b
        <Basepair(A103-A104)>

        # loading test basepair
        >>> b = instances.basepair()
        >>> print b
        <Basepair(A13-A12)>

        # primary axis of orientation, used to align to this basepair
        >>> print b.r()
        [[  1.00000001e+00   1.00000001e-04  -9.99999990e-09]
         [ -1.00000002e-04   1.00000001e+00  -9.99999990e-05]
         [ -1.99999995e-08   9.99999955e-05   1.00000000e+00]]

        # center of mass of basepair
        >>> print b.d()
        [ 0.1956032   0.69256601  0.0930465 ]

    """

    __slots__ = [
        "_res1_uuid",
        "_res2_uuid",
        "_r",
        "_d",
        "_sugars",
        "_name",
        "_bp_type",
        "_uuid"]

    def __init__(self, res1_uuid, res2_uuid, r, d, sugars, name,
                 bp_type=None, bp_uuid=None):
        self._res1_uuid, self._res2_uuid = res1_uuid, res2_uuid
        self._r = r
        self._d = d
        self._sugars = sugars
        self._name = name
        self._bp_type = bp_type
        self._uuid = bp_uuid

        if self._bp_type is None:
            self._bp_type = "c..."

        if self._uuid is None:
            self._uuid = uuid.uuid1()

    @classmethod
    def from_str(cls, s, res1_uuid, res2_uuid):
        spl = s.split(";")
        d = basic_io.str_to_point(spl[0])
        r = basic_io.str_to_matrix(spl[1])
        sugars = basic_io.str_to_points(spl[2])

        return cls(res1_uuid, res2_uuid, r, d, sugars, spl[3], spl[4])

    @classmethod
    def copy(cls, bp):
        sugars = [np.copy(bp._sugars[0]), np.copy(bp._sugars[1])]
        return cls(bp._res1_uuid, bp._res2_uuid, np.copy(bp._r), np.copy(bp._d),
                   sugars, bp._name, bp._bp_type, bp._uuid)

    @classmethod
    def copy_with_new_uuids(cls, bp, res1_uuid, res2_uuid):
        sugars = [np.copy(bp._sugars[0]), np.copy(bp._sugars[1])]
        return cls(res1_uuid, res2_uuid, np.copy(bp._r), np.copy(bp._d),
                   sugars, bp._name, bp._bp_type, bp_uuid=uuid.uuid1())

    def __repr__(self):
          return "<Basepair("+self._name + ")>"

    def partner(self, r_uuid):
        """
        get the other basepairing partner of a residue will throw an error
        if the supplied residue is not contained within this basepair

        :param res: the residue that you want to get the partner of
        :type res: secondary_structure.Residue object
        """

        if   r_uuid == self._res1_uuid:
            return self.res2_uuid
        elif r_uuid == self._res2_uuid:
            return self.res1_uuid
        else:
            raise exceptions.BasepairException(
                "call partner with a residue not in basepair")

    def get_transforming_r_and_t(self, r, t, sugars):
        """
        get a rotation matrix and translation that describes the tranformation
        between the rotation, translation to THIS BasepairState.

        :param r: Another orientation matrix from another basepair
        :param t: The origin of another basepair
        :param sugars: the c1' coords of another basepair

        :type r: np.array
        :type t: np.array
        :type sugars: list of two np.arrays

        :return: rotation and translation that defines transformation betwen
            both states
        """

        r1 = self.r
        r2 = r
        r_trans = util.unitarize(r1.T.dot(r2))
        t_trans = -t

        new_sugars_2 = np.dot(sugars, r_trans.T) + t_trans + self.d

        if sugars is not None:
            diff = (((self.sugars[0] - new_sugars_2[0]) +
                     (self.sugars[1] - new_sugars_2[1]))/2)
        else:
            diff = 0
        return r_trans, t_trans+diff

    def get_transforming_r_and_t_w_state(self, state):
        """
        wrapper for get_transforming_r_and_t using another basepair state
        instead of specifying each component explicitly.

        :param state: the basepair state you would like to get a transformation
            to align to the basepair state calling this function
        :type state: BasepairState

        """

        return self.get_transforming_r_and_t(state._r, state._d, state._sugars)

    def get_transformed_state(self, r, t):
        """
        get new orientation, origin and sugar coordinates after transforming
        with suplied rotation and translation.

        :param r: supplied rotation matrix
        :param t: supplied translation

        :type r: np.array
        :type t: np.array

        :return: new orientation, origin and sugar coorindates of this basepair
            state after transformation

        """

        r_T = r.T

        new_r = util.unitarize(np.dot(self._r, r_T))
        new_sugars = np.dot(self._sugars, r_T) + t
        new_origin = np.dot(self._d, r_T) + t

        return new_r, new_origin, new_sugars

    def transform(self, t):
        r_T = t.rotation().T

        new_r = util.unitarize(np.dot(self._r, r_T))
        new_sugars = []
        for s in self._sugars:
            new_sugars.append(np.dot(s, r_T) + t.translation())
        new_origin = np.dot(self._d, r_T) + t.translation()

        self._r = new_r
        self._d = new_origin
        self._sugars = new_sugars

    def fast_transform(self, r, t):
        r_T = r

        new_r = util.unitarize(np.dot(self._r, r_T))
        new_sugars = []
        for s in self._sugars:
            new_sugars.append(np.dot(s, r_T) + t)
        new_origin = np.dot(self._d, r_T) + t

        self._r = new_r
        self._d = new_origin
        self._sugars = new_sugars

    def move(self, p):
        self._sugars[0] += p
        self._sugars[1] += p
        self._d += p

    def diff(self, bp):
        diff = util.distance(self.d, bp.d)
        diff += self._rot_diff(bp) * 2
        return diff

    def _rot_diff(self, bp):
        r_diff = util.matrix_distance(self.r, bp.r)
        bp.flip()
        r_diff_2 = util.matrix_distance(self.r, bp.r)
        bp.flip()
        if r_diff > r_diff_2:
            r_diff = r_diff_2
        return r_diff

    def to_str(self):
        s  = basic_io.point_to_str(self._d) + ";"
        s += basic_io.matrix_to_str(self._r) + ";"
        s += basic_io.points_to_str(self._sugars) + ";"
        s += self._name + ";" + self._bp_type + ";"
        return s

    def flip(self):
        self.r[1] = -self.r[1]
        self.r[2] = -self.r[2]

    @property
    def r(self):
        return self._r

    @property
    def d(self):
        return self._d

    @property
    def sugars(self):
        return self._sugars

    @property
    def res1_sugar(self):
        return self._sugars[0]

    @property
    def res2_sugar(self):
        return self._sugars[1]

    @property
    def bp_type(self):
        return self._bp_type

    @property
    def uuid(self):
        return self._uuid

    @property
    def res1_uuid(self):
        return self._res1_uuid

    @property
    def res2_uuid(self):
        return self._res2_uuid

    @property
    def name(self):
        return self._name


class Motif(primitives.rna_structure.RNAStructure):
    __slots__ = [
        "_structure",
        "_basepairs",
        "_ends",
        "_name",
        "_score",
        "_end_ids",
        "_block_end_add",
        "_uuid",
        "_mtype",
        "_score"
    ]

    def __init__(self, structure, basepairs, ends, end_ids, name, mtype, score,
                 block_end_add=0, m_uuid=None):

        super(self.__class__, self).__init__(structure, basepairs, ends,
                                             end_ids, name)

        self._block_end_add = block_end_add
        self._mtype = mtype
        self._score = score
        self._uuid = m_uuid

        if self._uuid is None:
            self._uuid = uuid.uuid1()

    @classmethod
    def from_str(cls, s):
        spl = s.split("&")
        name = spl[0]
        score = float(spl[1])
        block_end_add = int(spl[2])
        mtype = int(spl[3])
        struc = Structure.from_str(spl[4])
        bp_strs = spl[5].split("@")
        bps = []
        for bp_str in bp_strs[:-1]:
            bps.append(bp_from_str(struc, bp_str))
        ends = [ bps[int(i)] for i in spl[6].split() ]
        end_ids = spl[7].split()
        return cls(struc, bps, ends, end_ids, name, mtype, score,
                   block_end_add)

    @classmethod
    def copy(cls, m, new_uuid=0):
        s = Structure.copy(m._structure, new_uuid)
        basepairs = []
        ends = []

        for bp in m._basepairs:
            if new_uuid:
                bp_res = m.get_bp_res(bp)
                res1 = s.get_residue(num=bp_res[0].num, chain_id=bp_res[0].chain_id,
                                     i_code=bp_res[0].i_code)
                res2 = s.get_residue(num=bp_res[1].num, chain_id=bp_res[1].chain_id,
                                     i_code=bp_res[1].i_code)
                bp = Basepair.copy_with_new_uuids(bp, res1.uuid, res2.uuid)
                basepairs.append(bp)
            else:
                basepairs.append(Basepair.copy(bp))

        for end in m._ends:
            i = m._basepairs.index(end)
            ends.append(basepairs[i])

        return cls(s, basepairs, ends, m._end_ids, m._name, m._mtype, m._score,
                   m._block_end_add, m._uuid)

    def to_str(self):
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
        return s

    def move(self, p):
        self._structure.move(p)
        for bp in self._basepairs:
            bp.move(p)

    def transform(self, t):
        self._structure.transform(t)
        for bp in self._basepairs:
            bp.transform(t)

    def fast_transform(self, t):
        r = t.rotation().T
        trans = t.translation()
        self._structure.fast_transform(r, trans)
        for bp in self._basepairs:
            bp.fast_transform(r, trans)

    @property
    def mtype(self):
        return self._mtype

    @property
    def uuid(self):
        return self._uuid

    @property
    def block_end_add(self):
        return self._block_end_add

    @property
    def score(self):
        return self._score

class MotifEnsemble(Ensemble):
    __slots__ = [
        "_end_id",
        "_members",
        "_block_end_add"
    ]

    def __init__(self, motifs, energies):
        if len(motifs) != len(energies):
            raise ValueError("must supply the same number of motifs and energies")

        members = []
        for i, m in enumerate(motifs):
            ms = EnsembleMember(m, energies[i])
            members.append(ms)

        members.sort(key=lambda x: x.energy, reverse=False)
        super(self.__class__, self).__init__(members)

    def copy(self):
        mes_copy = MotifStateEnsemble()

        mes_copy.id =self.id
        mes_copy.block_end_add = self.block_end_add
        members = []
        for mem in self.members:
            members.append(mem.copy())
        mes_copy.members = members

        return mes_copy

    def to_str(self):
        s = self.id + "{" + str(self.block_end_add) + "{"
        for ms in self.members:
            s += ms.to_str() + "{"
        return s


class MotifLibrary(object):
    def __init__(self, name, motifs=None):
        self._name = name
        self._motifs = motifs
        if self._motifs is None:
            self._motifs = []

    def __iter__(self):
        return self._motifs.__iter__()

    @property
    def name(self):
        return self._name

    @property
    def motifs(self):
        return self._motifs



def align_motif_state(ref_bp_state, ms):
    r, t = ref_bp_state.get_transforming_r_and_t_w_state(ms.get_end(0))
    t += ref_bp_state.d

    trans  = transform.Transform(r, t)

    #ms.transform(trans)
    ms.fast_transform(trans)

def get_aligned_motif_state(ref_bp_state, ms, new_uuid=1):
    r, t = ref_bp_state.get_transforming_r_and_t_w_state(ms.get_end(0))
    t += ref_bp_state.d

    ms_copy = Motif.copy(ms, new_uuid=new_uuid)
    trans  = transform.Transform(r, t)

    ms_copy.transform(trans)

    return ms_copy


def bp_from_str(struc, s):
    bp_spl = s.split(";")
    r2_info = bp_spl.pop().split("|")
    r1_info = bp_spl.pop().split("|")
    bp_str = ";".join(bp_spl)
    res1 = struc.get_residue(int(r1_info[0]), r1_info[1], r1_info[2])
    res2 = struc.get_residue(int(r2_info[0]), r2_info[1], r2_info[2])
    return Basepair.from_str(bp_str, res1.uuid, res2.uuid)


















