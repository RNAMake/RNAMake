import util
import basic_io
import settings
import motif
import numpy as np
import uuid
import exceptions


class Basepair(object):
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

    __slots__ = ["res1", "res2", "r", "bp_type", "atoms", "bp_state", "uuid"]

    def __init__(self, res1, res2, r, bp_type="c..."):
        self.res1, self.res2 = res1, res2
        self.atoms = self._get_atoms()
        self.bp_state = self._get_state(r)
        self.bp_type = bp_type
        self.uuid = uuid.uuid1()

    def __repr__(self):
          return "<Basepair("+self.name()+ ")>"

    def _get_atoms(self):
        """
        helper function to collect all non null atoms from both resiudes in
        basepair

        :return: all atoms that belong to basepair
        :rtype: list of atom.Atoms
        """

        atoms = []
        for a in self.res1.atoms:
            if a is not None:
                atoms.append(a)
        for a in self.res2.atoms:
            if a is not None:
                atoms.append(a)
        return atoms

    def get_base_atoms(self):
        atoms = []
        for a in self.res1.atoms[12:]:
            if a is not None:
                atoms.append(a)
        for a in self.res2.atoms[12:]:
            if a is not None:
                atoms.append(a)
        return atoms

    def _get_state(self, r):
        """
        helper function to setup basepair state object that keeps track of the
        rotation, center and c1' sugar atoms used to align to this basepair

        :param r: orientation matrix defining principle axes of coordinate
            system
        :type r: np.array

        :return: basepair state object for this basepair
        :rtype: BasepairState
        """

        d = util.center(self.atoms)
        sugars = [self.res1.get_atom("C1'").coords,
                  self.res2.get_atom("C1'").coords]
        return BasepairState(r, d, sugars)

    def residues(self):
        """
        returns both residues for easy iteration

        :returns: both residues in base
        :rtype: list of residue.Residues
        """
        return [self.res1, self.res2]

    def c1_prime_coords(self):
        """
        gets the c1' atom coordinates for both residues. These coordinates are
        used to fine tune the alignment between basepairs

        :return: list of c1' coords
        :rtype: list of np.arrays
        """

        return [self.res1.get_atom("C1'").coords,
                self.res2.get_atom("C1'").coords]

    def partner(self, res):
        """
        get the other basepairing partner of a residue will throw an error
        if the supplied residue is not contained within this basepair

        :param res: the residue that you want to get the partner of
        :type res: Residue object
        """

        if res == self.res1:
            return self.res2
        elif res == self.res2:
            return self.res1
        else:
            raise exceptions.BasepairException(
                "called get_partner with residue not in this not in this "
                "basepair")

    def name(self):
        """
        get name of basepair: which is the combined name of both residues
        seperated by a "-". The residue with the lower res number should
        come first

        :return: name of basepair
        :rtype: str

        :examples:

        ..  code-block:: python

            # build basepair from stratch
            >>> from rnamake.unittests import instances
            >>> b = instances.basepair()
            >>> print b.res1
            <Residue('G13 chain A')>

            >>> print b.res2
            <Residue('C12 chain A')>

            >>> print b.name()
            A12-A13
        """

        res1_name = self.res1.chain_id+str(self.res1.num)+str(self.res1.i_code)
        res2_name = self.res2.chain_id+str(self.res2.num)+str(self.res2.i_code)

        if self.res1.chain_id < self.res2.chain_id:
            return res1_name+"-"+res2_name
        if self.res1.chain_id > self.res2.chain_id:
            return res2_name+"-"+res1_name

        if self.res1.num < self.res2.num:
            return res1_name+"-"+res2_name
        else:
            return res2_name+"-"+res1_name

    def flip(self):
        """
        wrapper for BasepairState flip. There is probably no reasons to ever
        call this. See the implementation in BasepairState

        :return: None
        """

        self.bp_state.flip()

    def copy(self):
        """
        creates a deep copy of this basepair object
        """
        cbp = Basepair(self.res1, self.res2, self.bp_state.r)
        cbp.bp_state = self.bp_state.copy()
        cbp.bp_state.sugars = [self.res1.get_atom("C1'").coords,
                               self.res2.get_atom("C1'").coords]
        cbp.bp_type = self.bp_type
        cbp.uuid = self.uuid
        return cbp

    def state(self):
        """
        gets the state of this basepair and makes sure that the
        center and sugar positions are updated
        """
        d = util.center(self.atoms)
        self.bp_state.d = d
        self.bp_state.sugars = self.c1_prime_coords()
        return self.bp_state

    def r(self):
        """
        gets the orientation matrix of this basepair which is stored in the
        the BasepairState instance.

        :return: orientation matrix of basepair
        :rtype: np.array
        """

        return self.bp_state.r

    def d(self):
        """
        gets the center of the atoms in this basepair.

        :return: center of mass of basepair
        :rtype: np.array
        """

        return util.center(self.atoms)

    def to_str(self):
        """
        stringify basepair object
        """
        str1 = self.res1.chain_id+str(self.res1.num)+str(self.res1.i_code)
        str2 = self.res2.chain_id+str(self.res2.num)+str(self.res2.i_code)

        s = str1+"-"+str2 + "," + self.bp_state.to_str() + "," + \
            self.bp_type
        return s

    def to_pdb_str(self, acount=1, return_acount=0):
        """
        creates a PDB string formatted verision of this Basepair object.

        :param acount: current atom index, default: 1
        :param return_acount: final atom index after current atoms, default: 0

        :type  acount: int
        :type  return_acount: int

        :return: str
        """
        s = ""
        for r in self.residues():
            r_str, acount = r.to_pdb_str(acount, 1)
            s += r_str
        if return_acount:
            return s, acount
        else:
            return s

    def to_pdb(self, fname="basepair.pdb"):
        """
        write basepair object to pdb file

        :param fname" the file you want to write the basepair too
        :type fname: str
        """
        f = open(fname, "w")
        f.write ( self.to_pdb_str() )
        f.close()

    def transform(self, t):

        r_T = t.rotation().T
        for a in self.atoms:
            a.coords = np.dot(a.coords, r_T) + t.translation()

        transformed = np.dot(self.state().r, r_T)
        self.state().r = transformed

    def move(self, p):
        for a in self.atoms:
            a.coords += p

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


    def __repr__(self):
        return """\n\
        <BasepairState>\n\
        ...Frame...\n\
        ||\t%f\t%f\t%f\t%f\t||\n\
        ||\t%f\t%f\t%f\t%f\t||\n\
        ||\t%f\t%f\t%f\t%f\t||\n\
        ||\t%f\t%f\t%f\t%f\t||\n\
        ...Sugars...\n\
        [%f\t%f\t%f\t]\n\
        [%f\t%f\t%f\t]\n\
        """%\
       (self.r[0][0],self.r[1][0],self.r[2][0],self.d[0],
        self.r[0][1], self.r[1][1], self.r[2][1], self.d[1],
        self.r[0][2], self.r[1][2], self.r[2][2], self.d[2],
        0,0,0,1,
        self.sugars[0][0],self.sugars[0][1],self.sugars[0][2],
        self.sugars[1][0], self.sugars[1][1], self.sugars[1][2]
        )

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

