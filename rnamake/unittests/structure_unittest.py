import unittest
import os
import rnamake.structure
import rnamake.settings
import util

class StructureUnittest(unittest.TestCase):

    def test_build_chains(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6.pdb"
        struct = util.supress_log_output(rnamake.structure.Structure,
                                         path)

        if len(struct.chains) != 1:
            self.fail("did not get the expected number of chains")

    def _generate_build_chains_test(self, f_path):
        try:
            import redesign.motif
        except:
            self.skipTest("cannot import old resdesign package")

        path = "/Users/josephyesselman/projects/REDESIGN/resources/non-redundant-rnas"
        if not os.path.isdir(path):
            self.skipTest("non-redunant-rna path is not set cannot do this test")

        dirs = []
        for x in os.listdir(path):
            if os.path.isdir(path + "/" + x):
                dirs.append(x)

        f = open(f_path,"w")

        for d in dirs:
            try:
                m = redesign.motif.Motif(path + "/" + d)
            except:
                continue
            f.write(d + ",")
            for c in m.structure.chains:
                start = c.residues[0]
                end = c.residues[-1]

                start_key = start.name + " " + str(start.num) + " " + start.chain_id
                end_key = end.name + " " + str(end.num) + " " + end.chain_id
                f.write(start_key + " " + end_key + ",")
            f.write("\n")
        f.close()

    def test_build_chains_all(self):
        if util.UnittestState == util.UnittestType.BASIC:
            self.skipTest("Unittest State set to BASIC")

        path = rnamake.settings.UNITTEST_PATH + "resources/build_chains.dat"
        if not os.path.isfile(path):
            self._generate_build_chains_test(path)

        f = open(path)
        lines = f.readlines()
        f.close()

        skip_structs = "3SLQ".split()

        path = "/Users/josephyesselman/projects/REDESIGN/resources/non-redundant-rnas"
        for l in lines:
            spl = l.split(",")
            print spl[0]
            if spl[0] in skip_structs:
                continue
            pdb_path = path + "/" + spl[0] + "/" + spl[0] + ".pdb"
            struct = util.supress_log_output(rnamake.structure.Structure,
                                             pdb_path)
            for ckey in spl[1:-1]:
                found = 0
                for c in struct.chains:
                    start = c.first()
                    end = c.last()
                    start_key = start.name + " " + str(start.num) + " " + start.chain_id
                    end_key = end.name + " " + str(end.num) + " " + end.chain_id
                    key = start_key + " " + end_key
                    if key == ckey:
                        found = 1
                        break
                if not found:
                    print pdb_path,ckey
                    self.fail()

    def test_get_residue(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6.pdb"
        struct = util.supress_log_output(rnamake.structure.Structure,
                                         path)
        try:
            struct.get_residue()
            self.fail()
        except ValueError:
            pass
        except:
            self.fail("got an error I did not expect")

        res = struct.get_residue(num=107)
        if res is None:
            self.fail("should not of gotten an error")

        res = struct.get_residue(num=1000)
        if res is not None:
            self.fail("should not of gotten an error")

    def test_residues(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6.pdb"
        struct = util.supress_log_output(rnamake.structure.Structure,
                                         path)

        residues = struct.residues()
        if len(residues) != 157:
            self.fail()

    def test_atoms(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6.pdb"
        struct = util.supress_log_output(rnamake.structure.Structure,
                                         path)
        atoms = struct.atoms()
        if len(atoms) != 3357:
            self.fail()

def main():
    unittest.main()

if __name__ == '__main__':
    main()
