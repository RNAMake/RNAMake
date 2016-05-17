import unittest
import sys
import os
import rnamake.pdb_parser
import rnamake.residue_type
import rnamake.residue
import rnamake.settings
import rnamake.atom
import util
import warnings
from rnamake.util import distance
from rnamake import exceptions


class PdbParserUnittest(unittest.TestCase):

    def test_parse(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6.pdb"
        residues = rnamake.pdb_parser.parse(path)

    def test_warnings(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/pdbs/p4p6_error_1.pdb"

        #catch warning for missing atom name
        message = "line 2: no atomname detected"
        with warnings.catch_warnings(record=True) as w:
            rnamake.pdb_parser.parse(path)
            self.failUnless(str(w[-1].message) == message)

    def _get_residues_from_prody(self, prody_structure):
        """
        taken from old REDESIGN code, coverts a prody structure into, new
        residue objects. returns a list of residue objects
        """
        hv = prody_structure.getHierView()
        p_chains = list(hv)
        residues = []
        for pc in p_chains:
            for pr in list(pc):
                rtype = rnamake.residue_type.get_rtype(pr.getResname())
                if rtype is None:
                    continue
                res = rnamake.residue.Residue(rtype, pr.getResname(),
                                              pr.getResnum(), pr.getChid(),
                                              pr.getIcode())
                atoms = []
                for pa in list(pr):
                    atoms.append(rnamake.atom.Atom(pa.getName(), pa.getCoords()))
                res.setup_atoms(atoms)
                residues.append(res)

        return residues

    def _are_residues_the_same(self, residues1, residues2):
        r_hash_1 = { repr(r) : r for r in residues1 }
        r_hash_2 = { repr(r) : r for r in residues2 }

        for k,res1 in r_hash_1.iteritems():
            res2 = r_hash_2.get(k)
            if res2 is None:
                print k
                self.fail()
            for i in range(len(res1.atoms)):
                if   res1.atoms[i] is None and res2.atoms[i] is None:
                    continue
                elif res1.atoms[i] is None and res2.atoms[i] is not None:
                    self.fail("in res: " + res1 + "atom :" + \
                              res1.rtype.atom_map[i] + " is None in " + \
                              "res1 but not res2" )
                elif res1.atoms[i] is not None and res2.atoms[i] is None:
                     self.fail("in res: " + res1 + "atom :" + \
                              res1.atoms[i].name + " is None in " + \
                              "res2 but not res1" )
                d = distance(res1.atoms[i].coords,res2.atoms[i].coords)
                if d > 0.01:
                    print res1.atoms[i].name, res1
                    print res1.atoms[i].coords
                    print res2.atoms[i].coords
                    self.fail("distance between same atoms too great")

    def test_parse_compare(self):
        try:
            import prody
        except:
            self.skipTest("cannot import prody")
        prody.LOGGER.verbosity = 'none'
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6.pdb"
        new_residues = util.supress_log_output(rnamake.pdb_parser.parse,path)

        prody_structure = prody.parsePDB(path)
        old_residues = self._get_residues_from_prody(prody_structure)

        self._are_residues_the_same(new_residues, old_residues)

    def _test_parse_compare_all(self):
        #this test can take a while, default set to not run
        #if util.UnittestState == util.UnittestType.BASIC:
        #    self.skipTest("test_parse_compare_all is not a basic test")

        try:
            import prody
        except:
            self.skipTest("cannot import prody")

        path = "/Users/josephyesselman/projects/REDESIGN/resources/non-redundant-rnas"
        if not os.path.isdir(path):
            self.skipTest("non-redunant-rna path is not set cannot do this test")

        #2ZUF: random unconnected C that prody doesnt process
        #3OK4: formatting is really weird rexamine
        skip_dirs = "2ZUF 3OK4".split()

        dirs = []
        for x in os.listdir(path):
            if os.path.isdir(path + "/" + x):
                dirs.append(x)

        for d in dirs:
            if d in skip_dirs:
                continue

            pdb_path = path + "/" + d + "/" + d + ".pdb"
            print pdb_path
            try:
                with warnings.catch_warnings(record=True) as w:
                    new_residues = rnamake.pdb_parser.parse(pdb_path)
            except exceptions.PDBParserException:
                continue
            prody_structure = prody.parsePDB(pdb_path)
            old_residues = self._get_residues_from_prody(prody_structure)

            self._are_residues_the_same(new_residues, old_residues)


def main():
    unittest.main()

if __name__ == '__main__':
    #if len(sys.argv) > 1:
    #    util.UnittestState = int(sys.argv[1])
    #    sys.argv.pop()

    main()
