import unittest
import build
import numpy as np
import rnamake.sqlite_library as sqlite_library
import rnamake.secondary_structure
import rnamake.motif_factory as motif_factory
import rnamake.motif_tree as motif_tree
import rnamake.ss_tree as ss_tree
import rnamake.resource_manager as rm
import rnamake.util as util


class SqliteLibraryUnittest(unittest.TestCase):

    def test_creation(self):
        mlib = sqlite_library.MotifSqliteLibrary("twoway")

    def test_get(self):
        mlib = sqlite_library.MotifSqliteLibrary("ideal_helices")
        m1 = mlib.get(name="HELIX.IDEAL.6")
        m2 = mlib.get(name='HELIX.IDEAL.6', end_name='A8-A9')
        m3 = mlib.get(end_id='GGGGGGGG_LLLLLLLL_CCCCCCCC_RRRRRRRR')

        if m1 is None or m2 is None or m3 is None:
            self.fail("something wrong in get()")

        mlib = sqlite_library.MotifStateSqliteLibrary('ideal_helices')
        ms1 = mlib.get(name="HELIX.IDEAL.6")
        ms2 = mlib.get(name='HELIX.IDEAL.6', end_name='A8-A9')
        ms3 = mlib.get(end_id='GGGGGGGG_LLLLLLLL_CCCCCCCC_RRRRRRRR')

        if ms1 is None or ms2 is None or ms3 is None:
            self.fail("something wrong in get()")

        me_lib = sqlite_library.MotifEnsembleSqliteLibrary('bp_steps')
        me = me_lib.get(name='GG_LL_CC_RR')

        #print len(me.members)

    def test_load_all(self):
        mlib = sqlite_library.MotifSqliteLibrary("ideal_helices")
        mlib.load_all()
        if  len(mlib.data) == 0:
            self.fail("something wrong with load_all()")

    def test_get_multi(self):
        mlib = sqlite_library.MotifSqliteLibrary("twoway")
        m = mlib.get_random()
        motifs1 = mlib.get_multi(name=m.name)
        motifs2 = mlib.get_multi(end_id=m.end_ids[0])

    def test_bp_steps(self):
        mlib = sqlite_library.MotifSqliteLibrary("bp_steps")
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

            mt.remove_node()
            mt.remove_node()

        builder = build.BuildMotifTree(libs=[mlib])

        """for i in range(10):
            mt = builder.build(10)
            if len(mt) != 10:
                for n in mt:
                    print n.data.name
                print len(mt)
                mt.write_pdbs()
                self.fail("random generation of helices failed")"""

    def test_twoways(self):

        builder = build.BuildMotifTree()

        for i in range(10):
            mt = builder.build(3)
            if len(mt) != 3:
                for n in mt:
                    print n.data.name
                print len(mt)
                mt.write_pdbs()
                self.fail("random generation of twoways failed")

    def _get_random_motif(self, size=5):
        mlib = sqlite_library.MotifSSIDSqliteLibrary("twoway")
        m = mlib.get_random()
        while len(m.residues()) > 5:
            m = mlib.get_random()

        return m

    def test_get_best(self):
        mlib = sqlite_library.MotifSqliteLibrary("twoway")
        #ss_id1 = "AC_LL_GGGU_RUUR"
        ss_id1 = "CGG_LLL_CAG_RRR"

        matches = mlib.get_best_matches(ss_id1)
        #print len(matches)


    def _test_end_id(self):
        mlib = rm.manager.mlibs['tcontact']
        mlib.load_all()
        keep = "103 109 111 117 120".split()
        for i, m in enumerate(mlib.all()):
            if len(m.ends) != 3:
                continue
            if str(i) in keep:
                print m.name

    def _test_specific(self):
        ss_tree_1 = rnamake.secondary_structure.ss_id_to_ss_tree("AC_LL_GGU_RUR")
        ss_tree_2 = rnamake.secondary_structure.ss_id_to_ss_tree("AGAC_LUUL_GGGU_RUUR")

        for i, n in enumerate(ss_tree_1):
            print n.data.seq(), ss_tree_2.get_node(i).data.seq()

        print ss_tree.compare_ss_tree(ss_tree_1, ss_tree_2)

    def _test_get_by_topology(self):
        mlib = sqlite_library.MotifSSIDSqliteLibrary("twoway")
        motifs = mlib.get_by_topology([1,0])

    def test_end_orientation(self):
        mlib = sqlite_library.MotifSqliteLibrary("bp_steps")
        mlib.load_all()

        for i,m in enumerate(mlib.all()):
            res = m.residues()
            vec1 = res[0].get_atom("C2'").coords - res[0].get_atom("O4'").coords
            vec2 = res[-1].get_atom("C2'").coords - res[-1].get_atom("O4'").coords
            result = np.dot(vec1, vec2)
            if result > 2:
                print i, result
                m.to_pdb("motif."+str(i)+".pdb")

    def _test_bp_step_ensembles(self):
        me_lib = sqlite_library.MotifEnsembleSqliteLibrary("bp_steps")
        me = me_lib.get(name="GG_LL_CC_RR")

        for i, mem in enumerate(me.members):
            mem.motif.to_pdb("motif."+str(i)+".pdb")
            print i, mem.energy

    def _test_bp_steps_2(self):
        mlib = sqlite_library.MotifSqliteLibrary("bp_steps")
        mlib.load_all()
        for m in mlib.all():
            print m.name

    def test_get_1_0(self):
        mlib = sqlite_library.MotifSqliteLibrary("twoway")
        mlib.load_all()

        motifs = []
        end_name = ""
        for i,m in enumerate(mlib.all()):
            if len(m.residues()) != 5:
                continue
            if len(m.chains()[0].residues) != 3:
                continue
            if m.chains()[0].residues[1].name != 'A':
                continue

            if m.name == "TWOWAY.1DUQ.7":
                end_name = m.ends[0].name()
            motifs.append(m)
            #m.to_pdb("m."+str(i)+".pdb")
        print end_name


def main():
    unittest.main()

if __name__ == '__main__':
    main()
