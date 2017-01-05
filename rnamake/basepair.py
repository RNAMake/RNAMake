import util
import basic_io
import numpy as np
import uuid

import primitives

class Basepair(primitives.Basepair):
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
    def copy(cls, bp):
        sugars = [np.copy(bp._sugars[0]), np.copy(bp._sugars[1])]
        return cls(bp._res1_uuid, bp._res2_uuid, np.copy(bp._r), np.copy(bp._d),
                   sugars, bp._name, bp._bp_type, bp._uuid)

    def __repr__(self):
          return "<Basepair("+self.res1.chain_id+str(self.res1.num)+str(self.res1.i_code) +\
            "-" + self.res2.chain_id+str(self.res2.num)+str(self.res2.i_code) + ")>"

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

    def move(self, p):
        self._sugars[0] += p
        self._sugars[1] += p
        self._d += p

    def diff(self, bp):
        diff = util.distance(self.d(), bp.d())
        diff += self._rot_diff(bp) * 2
        return diff

    def _rot_diff(self, bp):
        r_diff = util.matrix_distance(self.r(), bp.r())
        bp.flip()
        r_diff_2 = util.matrix_distance(self.r(), bp.r())
        bp.flip()
        if r_diff > r_diff_2:
            r_diff = r_diff_2
        return r_diff

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

class BasepairState(object):
    """
    A small container class to hold the "State" of a basepair for finding
    matches in the database. The critical features are:

    :params r: Reference Frame of basepair
    :type r: np.array
    :params d: Center of mass of basepair
    :type d: np.array
    :params sugars: C1` atom coords for both residues in basepair
    :type sugars: List of Np.Arrays

    :attributes:
    `r` : np.Matrix
        Reference Frame of basepair
    `d` : Np.Array
        Center of mass of basepair
    `sugars` : List of Np.Arrays
        C1` atom coords for both residues in basepair

    :examples:

    ..  code-block:: python

        >>> from rnamake.unittests import instances
        >>> bp_state_1 = instances.basepairstate_random()
        >>> bp_state_2 = instances.basepairstate_random()

        # get rotation and translation defining the transformation between the
        # two states
        >>> r, t = bp_state_1.get_transforming_r_and_t_w_state(bp_state_2)
        >>> t += bp_state_1.d

        # align end2 to end1, notice this is two steps, this is done so
        # r and t can be applied to many states at the same time so they can
        # be all aligned together
        >>> new_r, new_d, new_sug = bp_state_2.get_transformed_state(r, t)
        >>> bp_state_2.set(new_r, new_d, new_sug)

        # both the orientation matrix and origins are nearly indentical now
        >>> print bp_state_1.r
        [[-0.43155336  0.24291737  0.86876513]
         [-0.61066224 -0.78751428 -0.08314381]
         [ 0.66396787 -0.56640305  0.4881949 ]]
        >>> print bp_state_2.r
        [[-0.43155336  0.24291737  0.86876513]
         [-0.61066224 -0.78751428 -0.08314381]
         [ 0.66396787 -0.56640305  0.4881949 ]]

        >>> print bp_state_1.d
        [  5.45492979  45.17079568  43.80163295]
        >>> print bp_state_2.d
        [  5.45492979  45.17079568  43.80163295]

    """

    __slots__ = ["r", "d", "sugars"]

    def __init__(self, r, d, sugars):
        self.r, self.d, self.sugars = r, d, sugars

    def flip(self):
        self.r[1] = -self.r[1]
        self.r[2] = -self.r[2]

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

        return self.get_transforming_r_and_t(state.r, state.d, state.sugars)

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

        new_r = util.unitarize(np.dot(self.r, r_T))
        new_sugars = np.dot(self.sugars, r_T) + t
        new_origin = np.dot(self.d, r_T) + t

        return new_r, new_origin, new_sugars

    def set(self, r, d, sug):
        """
        sets a new orientation matrx, origin and new c1' sugar coordinates,
        overriding the existing values stored in this instance.

        :param r: new orientation matrix
        :param d: new origin point
        :param sug: new c1' sugar coordinates

        :type r: np.array
        :type d: np.array
        :type sug: list of 2 np.arrays
        """

        self.r = r
        self.d = d
        self.sugars = sug

    def copy(self):
        """
        returns a deep copy of this BasepairState object
        """
        return BasepairState(np.array(self.r, copy=True),
                             np.array(self.d, copy=True),
                             [coord[:] for coord in self.sugars])

    def diff(self, state):
        """
        calculate the difference between this basepair and another. Difference
        is defined currently as euclidean distance between their origins and
        two times the matrix difference between the orientations of both
        states.

        :param state: the other state to compare the difference too
        :type state: BasepairState

        :return: difference between the states
        :rtype: float

        """

        diff  = util.distance(self.d, state.d)
        diff += self._rot_diff(state) * 2
        #diff += self._sugar_diff(state)
        return diff

    def _rot_diff(self, state):
        r_diff = util.matrix_distance(self.r, state.r)
        state.flip()
        r_diff_2 = util.matrix_distance(self.r, state.r)
        state.flip()
        if r_diff > r_diff_2:
            r_diff = r_diff_2
        return r_diff

    def _sugar_diff(self, state):
        diff_1 = util.distance(self.sugars[0], state.sugars[0]) + \
                 util.distance(self.sugars[1], state.sugars[1])
        diff_2 = util.distance(self.sugars[1], state.sugars[0]) + \
                 util.distance(self.sugars[0], state.sugars[1])
        if diff_1 > diff_2:
            diff_1 = diff_2
        return diff_1

    def to_str(self):
        """
        converts basepairstate into a string
        """

        s = basic_io.point_to_str(self.d) + ";" + \
            basic_io.matrix_to_str(self.r) + ";" +\
            basic_io.matrix_to_str(self.sugars)
        return s


def str_to_basepairstate(s):
    """
    convert stringified basepair back to a basepair object see
    basepairstate.to_str()
    """
    spl = s.split(";")
    d = basic_io.str_to_point(spl[0])
    r = basic_io.str_to_matrix(spl[1])
    sugars = basic_io.str_to_matrix(spl[2])
    return BasepairState(r, d, sugars)

