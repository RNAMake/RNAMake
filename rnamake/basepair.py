import util
import numpy as np

class Basepair(object):

    def __init__(self, res1, res2, r):
        self.res1, self.res2 = res1, res2
        self.atoms = self._get_atoms()
        self.bp_state = self._get_state(r)
        self.bp_type = "c..."
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
        sugars = [self.res1.get_atom("C1'"), self.res2.get_atom("C1'")]
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
        if   res == self.res1:
            return self.res2
        elif res == self.res2:
            return self.res1
        else:
            raise ValueError("called get_partner with residue not in this" +
                             "not in this basepair")

    def name(self):
         return self.res1.chain_id+str(self.res1.num)+str(self.res1.i_code) +\
         "-"+ self.res2.chain_id+str(self.res2.num)+str(self.res2.i_code)

    def flip(self,flip):
        if self.flip == flip:
            return
        else:
            self.state.flip()
            self.flipped = flip

    def copy(self):
        cbp = Basepair(self.res1, self.res2, self.bp_state.r)
        cbp.bp_state = self.bp_state.copy()
        cbp.bp_type = self.bp_type
        cbp.flipped = self.flipped
        cbp.designable = self.designable
        return cbp

class BasepairState(object):
    """
    A small container class to hold the "State" of a basepair for finding
    matches in the database. The critical features are:

    :params r: Reference Frame of basepair
    :type r: Np.Matrix
    :params d: Center of mass of basepair
    :type d: Np.Array
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

    def get_transforming_r_and_t(self, r, t ,sugars):

        r1 = self.r
        r2 = r
        r_trans = r1.T.dot(r2)
        t_trans = -t

        new_sugars_2 = np.dot(sugars, r_trans.T) + t_trans + self.d

        if sugars != None:
            diff = (((self.sugars[0] - new_sugars_2[0]) + \
                     (self.sugars[1] - new_sugars_2[1]))/2)
        else:
            diff = 0

        return r_trans,t_trans+diff

    def get_transformed_state(self,r,t):

        r_T = r.T

        new_r = unitarize(np.dot(self.r, r_T))
        new_sugars = np.dot(self.sugars, r_T) + t
        new_origin = np.dot(self.d, r_T) + t

        return new_r,new_origin,new_sugars

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

