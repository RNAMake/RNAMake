import unittest
import numpy as np
import sqlite3
import os

from rnamake import exceptions, sqlite_library, motif

class SqliteLibraryUnittest(unittest.TestCase):

    def test_creation(self):
        sqlite_library.MotifSqliteLibrary("twoway")

    def test_get(self):
        mlib = sqlite_library.MotifSqliteLibrary("ideal_helices")
        m1 = mlib.get(name="HELIX.IDEAL.6")
        m2 = mlib.get(name='HELIX.IDEAL.6', end_name='A1-A16')
        m3 = mlib.get(end_id='CCCCCCCC_LLLLLLLL_GGGGGGGG_RRRRRRRR')

        if m1 is None or m2 is None or m3 is None:
            self.fail("something wrong in get()")

        with self.assertRaises(exceptions.SqliteLibraryException):
            mlib.get(name="fake")

    def test_load_all(self):
        for k,v in sqlite_library.MotifSqliteLibrary.get_libnames().iteritems():
            mlib = sqlite_library.MotifSqliteLibrary(k)
            mlib.load_all(limit=10)
            if  len(mlib.data) == 0:
                self.fail("something wrong with load_all()")

    def test_get_multi(self):
        mlib = sqlite_library.MotifSqliteLibrary("twoway")
        m = mlib.get_random()
        motifs1 = mlib.get_multi(name=m.name)
        motifs2 = mlib.get_multi(end_id=m.get_end_id(0))

    #TODO check all move to intergration
    def _test_end_orientation(self):
        mlib = sqlite_library.MotifSqliteLibrary("bp_steps")
        mlib.load_all(10)

        for i,m in enumerate(mlib.all()):
            res = m.residues()
            vec1 = res[0].get_atom("C2'").coords - res[0].get_atom("O4'").coords
            vec2 = res[-1].get_atom("C2'").coords - res[-1].get_atom("O4'").coords
            result = np.dot(vec1, vec2)
            if result > 2:
                self.fail("wrong end direction")

    def test_do_libs_exist(self):
        for k in sqlite_library.MotifSqliteLibrary.get_libnames().keys():
            try:
                sqlite_library.MotifSqliteLibrary(k)
            except exceptions.SqliteLibraryException:
                self.fail("could not load a standard motif library!")

        """for k in sqlite_library.MotifStateSqliteLibrary.get_libnames().keys():
            try:
                sqlite_library.MotifStateSqliteLibrary(k)
            except exceptions.SqliteLibraryException:
                self.fail("could not load a standard motif_state library!")

        for k in sqlite_library.MotifEnsembleSqliteLibrary.get_libnames().keys():
            try:
                sqlite_library.MotifEnsembleSqliteLibrary(k)
            except exceptions.SqliteLibraryException:
                self.fail("could not load a standard motif_ensemble library!")

        for k in sqlite_library.MotifStateEnsembleSqliteLibrary.get_libnames().keys():
            try:
                sqlite_library.MotifStateEnsembleSqliteLibrary(k)
            except exceptions.SqliteLibraryException:
                self.fail("could not load a standard motif_state_ensemble library!")"""

    def test_build_sqlite_library(self):

        data = [['the_word', 0], ['the', 1], ['hello', 2]]
        keys = ['word', 'id']

        sqlite_library.build_sqlite_library("test.db", data, keys, 'id')

        conn = sqlite3.connect("test.db")
        r = conn.execute("SELECT * from data_table WHERE word=\'the_word\'").fetchone()

        self.failUnless(r[0] == 'the_word', "should be able to retrieve data")
        self.failUnless(int(r[1]) == 0, "should be able to retrieve data")

        os.remove("test.db")

    def test_new_res_uuids(self):
        mlib = sqlite_library.MotifSqliteLibrary("ideal_helices")
        m1 = mlib.get(name="HELIX.IDEAL")
        m2 = mlib.get(name="HELIX.IDEAL")
        self.failIf(m1.uuid == m2.uuid)

    def _test_get_all_bulges(self):
        mlib = sqlite_library.MotifSqliteLibrary("twoway")
        mlib.load_all()

        seen = []
        for m in mlib.all():
            if len(m.residues()) != 5:
                continue
            if m.name in seen:
                continue
            seen.append(m.name)

    def _test_get_1_0(self):
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
        #print end_name

    def _test_new_bp_steps(self):
        mlib = sqlite_library.MotifSqliteLibrary("new_bp_steps")
        mlib.load_all()

        end_indexes = []
        for m in mlib.all():
            if m.end_ids[0] not in end_indexes:
                end_indexes.append(m.end_ids[0])
            else:
                print m.name, m.end_ids[0]

    def _test_new_bp_steps(self):
        #mlib = sqlite_library.MotifSqliteLibrary("new_bp_steps")
        #mlib.load_all()

        m = rm.manager.get_bp_step("GC_LL_GC_RR")
        #m2 = mlib.get(end_id="GC_LL_GC_RR", end_name=m.ends[1].name())


        mg = motif_graph.MotifGraph()
        mg.add_motif(m)
        h = rm.manager.get_motif(name="HELIX.IDEAL.6")
        mg.add_motif(h)


        #m2.ends[0].flip()

        m_aligned = motif_graph.flip_alignment(m, 1)
        #m2.ends[0].flip()

        m.to_pdb("test.pdb")
        m_aligned.to_pdb("test_2.pdb")
        exit()

        mg2 = motif_graph.MotifGraph()
        mg2.add_motif(m_aligned)
        mg2.add_motif(rm.manager.get_motif(name="HELIX.IDEAL.6"))
        #mg.write_pdbs("one")
        #mg2.write_pdbs("two")

    def _test_helix(self):
        ms_lib = sqlite_library.MotifStateSqliteLibrary("ideal_helices")
        ms_lib.load_all()

        #m = rm.manager.get_motif(name="HELIX.IDEAL.7", end_name="A9-A10")
        #print m
        #exit()
        for m in ms_lib.all():
            print m.name, m.end_names[0]



def main():
    unittest.main()

if __name__ == '__main__':
    main()
