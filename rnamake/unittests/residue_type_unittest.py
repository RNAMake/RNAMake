import unittest
import logging
import numpy as np

from rnamake import atom, residue_type

class ResidueTypeUnittest(unittest.TestCase):

    def setUp(self):
        self.rts = residue_type.ResidueTypeSet()

    def test_creation(self):
       pass

    def test_get_rtype(self):
        gua_rtype = self.rts.get_rtype_by_resname("GUA")
        if gua_rtype is None:
            self.fail("cannot get GUA restype in test_get_rtype")

        none_rtype = self.rts.get_rtype_by_resname("GUA1")
        if none_rtype is not None:
            self.fail("did not return None properly in test_get_rtype")

        altname_rtype = self.rts.get_rtype_by_resname("I")
        if altname_rtype is None:
            self.fail("cannot get altname restype in test_get_rtype")

    def test_get_correct_atom_name(self):
        a1 = atom.Atom("O1P", np.array([0, 1, 2]))
        gua = self.rts.get_rtype_by_resname("GUA")
        self.failIf(gua.get_correct_atom_name(a1) is None)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
