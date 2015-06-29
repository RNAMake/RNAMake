import unittest
import build
import rnamake.sqlite_library as sqlite_library
import rnamake.secondary_structure
import rnamake.motif_factory as motif_factory
import rnamake.motif_tree as motif_tree
import rnamake.ss_tree as ss_tree

class SqliteLibraryUnittest(unittest.TestCase):


    def test_bp_steps(self):
        mlib = sqlite_library.MotifSSIDSqliteLibrary("bp_steps")
        mlib.load_all()

        base = motif_factory.factory.base_motif
        mt = motif_tree.MotifTree()
        mt.add_motif(base)

        for i, m in enumerate(mlib.all()):
            index = mt.add_motif(m)
            if index == -1:
                self.fail("something wrong with bp_step directionality")

            index = mt.add_motif(base)
            if index == -1:
                self.fail("something wrong with bp_step directionality 2")

            mt.remove_node(2)
            mt.remove_node(1)

        builder = build.BuildMotifTree(libs=[mlib])

        for i in range(10):
            mt = builder.build(100)
            if len(mt) != 100:
                for n in mt:
                    print n.data.name
                print len(mt)
                mt.write_pdbs()
                self.fail("random generation of helices failed")

    def test_twoways(self):

        builder = build.BuildMotifTree()

        for i in range(100):
            mt = builder.build(3)
            if len(mt) != 3:
                for n in mt:
                    print n.data.name
                print len(mt)
                mt.write_pdbs()
                self.fail("random generation of twoways failed")


    def test_creation(self):
        mlib = sqlite_library.MotifSSIDSqliteLibrary("twoway")

    def test_get(self):
        mlib = sqlite_library.MotifSSIDSqliteLibrary("twoway")
        m = mlib.get_random()
        while len(m.residues()) > 5:
            m = mlib.get_random()

        ss_id = m.end_ids[0]
        m1 = mlib.get(ss_id)

    def _get_random_motif(self, size=5):
        mlib = sqlite_library.MotifSSIDSqliteLibrary("twoway")
        m = mlib.get_random()
        while len(m.residues()) > 5:
            m = mlib.get_random()

        return m

    def test_get_best(self):
        mlib = sqlite_library.MotifSSIDSqliteLibrary("twoway")
        ss_id1 = "AC_LL_GGGU_RUUR"

        m = mlib.get_best_match(ss_id1)
        m.to_pdb("test.pdb")

    def _test_specific(self):
        ss_tree_1 = rnamake.secondary_structure.ss_id_to_ss_tree("AC_LL_GGU_RUR")
        ss_tree_2 = rnamake.secondary_structure.ss_id_to_ss_tree("AGAC_LUUL_GGGU_RUUR")

        for i, n in enumerate(ss_tree_1):
            print n.data.seq(), ss_tree_2.get_node(i).data.seq()

        print ss_tree.compare_ss_tree(ss_tree_1, ss_tree_2)


    def test_get_by_topology(self):
        mlib = sqlite_library.MotifSSIDSqliteLibrary("twoway")
        motifs = mlib.get_by_topology([1,0])




def main():
    unittest.main()

if __name__ == '__main__':
    main()
