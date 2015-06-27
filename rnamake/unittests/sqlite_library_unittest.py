import unittest
import rnamake.sqlite_library
import rnamake.secondary_structure
import rnamake.ss_tree as ss_tree
class SqliteLibraryUnittest(unittest.TestCase):

    def test_creation(self):
        mlib = rnamake.sqlite_library.MotifSSIDSqliteLibrary("twoway")

    def test_get(self):
        mlib = rnamake.sqlite_library.MotifSSIDSqliteLibrary("twoway")
        m = mlib.get_random()
        while len(m.residues()) > 5:
            m = mlib.get_random()

        ss_id = m.end_ids[0]
        m1 = mlib.get(ss_id)

    def _get_random_motif(self, size=5):
        mlib = rnamake.sqlite_library.MotifSSIDSqliteLibrary("twoway")
        m = mlib.get_random()
        while len(m.residues()) > 5:
            m = mlib.get_random()

        return m

    def test_get_best(self):
        mlib = rnamake.sqlite_library.MotifSSIDSqliteLibrary("twoway")
        ss_id1 = "AC_LL_GGGU_RUUR"

        m = mlib.get_best_match(ss_id1)
        print m.end_ids[0]
        m.to_pdb("test.pdb")

    def _test_specific(self):
        ss_tree_1 = rnamake.secondary_structure.ss_id_to_ss_tree("AC_LL_GGU_RUR")
        ss_tree_2 = rnamake.secondary_structure.ss_id_to_ss_tree("AGAC_LUUL_GGGU_RUUR")

        for i, n in enumerate(ss_tree_1):
            print n.data.seq(), ss_tree_2.get_node(i).data.seq()

        print ss_tree.compare_ss_tree(ss_tree_1, ss_tree_2)


    def test_get_by_topology(self):
        mlib = rnamake.sqlite_library.MotifSSIDSqliteLibrary("twoway")
        motifs = mlib.get_by_topology([1,0])




def main():
    unittest.main()

if __name__ == '__main__':
    main()
