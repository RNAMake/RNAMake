from . import util
from . import basic_io
import numpy as np


class Basepair(object):

    def __init__(self, res1, res2, r, bp_type="c..."):
        self.res1, self.res2 = res1, res2
        self.atoms = self._get_atoms()
        self.bp_state = self._get_state(r)
        self.bp_type = bp_type
        self.flipped = 0
        self.designable = 0

    def _get_atoms(self):
        atoms = []
        for a in self.res1.atoms:
            if a is not None:
                atoms.append(a)
        for a in self.res2.atoms:
            if a is not None:
                atoms.append(a)
        return atoms

    def _get_state(self, r):
        d = util.center(self.atoms)
        sugars = [self.res1.get_atom("C1'").coords,
                  self.res2.get_atom("C1'").coords]
        return BasepairState(r, d, sugars)

    def residues(self):
        """
        returns both residues for easy iteration
        """
        return [self.res1, self.res2]

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
            raise ValueError("called get_partner with residue not in this" +
                             "not in this basepair")

    def name(self):
        return self.res1.chain_id+str(self.res1.num)+str(self.res1.i_code) +\
            "-" + self.res2.chain_id+str(self.res2.num)+str(self.res2.i_code)

    def flip(self, flip):
        if self.flip == flip:
            return
        else:
            self.state.flip()
            self.flipped = flip

    def copy(self):
        """
        creates a deep copy of this basepair object
        """
        cbp = Basepair(self.res1, self.res2, self.bp_state.r)
        cbp.bp_state = self.bp_state.copy()
        cbp.bp_state.sugars =  [self.res1.get_atom("C1'").coords,
                                self.res2.get_atom("C1'").coords]
        cbp.bp_type = self.bp_type
        cbp.flipped = self.flipped
        cbp.designable = self.designable
        return cbp

    def state(self):
        """
        gets the state of this basepair and makes sure that the
        center and sugar positions are updated
        """
        d = util.center(self.atoms)
        self.bp_state.d = d
        return self.bp_state

    def to_str(self):
        """
        stringify basepair object
        """
        s = self.name() + "," + self.bp_state.to_str() + "," + \
            self.bp_type + "," + str(self.designable) + "," + \
            str(self.flipped)
        return s

    def to_pdb_str(self, acount=1, return_acount=0):
        """
        writes basepair object to a pdb formmatted string

        :param acount: the atom index start (optional)
        :param return_acount: whether the final atom index should be returned
            (optional)

        :type acount: int
        :type return_acount: int
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
        f = open(fname)
        f.write ( self.to_pdb_str() )
        f.close()

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

    Attributes
    ----------
    `r` : Np.Matrix
        Reference Frame of basepair
    `d` : Np.Array
        Center of mass of basepair
    `sugars` : List of Np.Arrays
        C1` atom coords for both residues in basepair

    """

    __slots__ = ["r", "d", "sugars"]

    def __init__(self, r, d, sugars):
        self.r, self.d, self.sugars = r, d, sugars

    def flip(self):
        self.r[1] = -self.r[1]
        self.r[2] = -self.r[2]

    def get_transforming_r_and_t(self, r, t, sugars):

        r1 = self.r
        r2 = r
        r_trans = r1.T.dot(r2)
        t_trans = -t

        new_sugars_2 = np.dot(sugars, r_trans.T) + t_trans + self.d

        if sugars is not None:
            diff = (((self.sugars[0] - new_sugars_2[0]) +
                     (self.sugars[1] - new_sugars_2[1]))/2)
        else:
            diff = 0

        return r_trans, t_trans+diff

    def get_transformed_state(self, r, t):

        r_T = r.T

        new_r = unitarize(np.dot(self.r, r_T))
        new_sugars = np.dot(self.sugars, r_T) + t
        new_origin = np.dot(self.d, r_T) + t

        return new_r, new_origin, new_sugars

    def set(self, r, d, sug):
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
