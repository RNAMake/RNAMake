import unittest
import random
import rnamake.motif_library_sqlite
import rnamake.motif_type
import rnamake.util

class MotifLibrarySqliteUnittest(unittest.TestCase):

    def test_creation(self):
        mtype = rnamake.motif_type.HELIX
        mlib = rnamake.motif_library_sqlite.MotifLibrarySqlite(mtype)

    def test_get_motif(self):
        mtype = rnamake.motif_type.HELIX
        mlib = rnamake.motif_library_sqlite.MotifLibrarySqlite(mtype)
        m = mlib.get_motif("HELIX.IDEAL")

        try:
            mlib.get_motif("HELIX")
            self.fail()
        except ValueError:
            pass
        except:
            self.fail("did not get the error I expected")

    def test_load_all(self):
        mtype = rnamake.motif_type.HELIX
        mlib = rnamake.motif_library_sqlite.MotifLibrarySqlite(mtype)
        mlib.load_all(limit=10)
        #for mname, m in mlib.mdict.iteritems():
        #    print mname, m

    def _build_test_tree(self, h_lib, twoway):
        mt = rnamake.motif_tree.MotifTree()
        size = 10
        pos = 0
        count = 0
        i = 0
        while i < size:
            if i % 2 == 0:
                pos = 0
            else:
                pos = 1

            if pos == 0:
                num = random.randint(1, 12)
                node = mt.add_motif(h_lib.get_motif("HELIX.IDEAL."+str(num)), end_flip=0, end_index=1)
            else:
                m = random.choice(twoway.motifs())
                spl = m.name.split("-")
                ei = int(spl[1])
                node = mt.add_motif(m, end_index=m.end_to_add, end_flip=0)


            if node is not None:
                i += 1

            count += 1
            if count > 1000:
                break

        return mt

    def _test_compare_2(self):
        mtype = rnamake.motif_type.HELIX
        h_lib = rnamake.motif_library_sqlite.MotifLibrarySqlite(mtype)
        path = "../resources/motif_libraries/twoway_aligned.db"
        twoways_new = rnamake.motif_library_sqlite.MotifLibrarySqlite(libpath=path)
        twoways_new.load_all()

        for m in twoways_new.motifs():
            print m.name
            mt = rnamake.motif_tree.MotifTree(sterics=1)
            spl = m.name.split("-")
            ei = int(spl[1])
            for i in range(1,13):
                mt.add_motif(h_lib.get_motif("HELIX.IDEAL."+str(i)),  end_flip=0, end_index=1)
                mt.add_motif(m, end_flip=0)
                mt.add_motif(h_lib.get_motif("HELIX.IDEAL."+str(i)),  end_flip=0, end_index=1)
                if len(mt.nodes) != 4:
                    print i, m.name, len(mt.nodes)

                mt.remove_node_level()


    def test_compare(self):
        mtype = rnamake.motif_type.HELIX
        h_lib = rnamake.motif_library_sqlite.MotifLibrarySqlite(mtype)
        twoways = rnamake.motif_library.unique_twoway_lib()
        path = "../resources/motif_libraries/twoway_aligned.db"
        twoways_new = rnamake.motif_library_sqlite.MotifLibrarySqlite(libpath=path)
        twoways_new.load_all()

        for i in range(100):
            mt2 = rnamake.motif_tree.MotifTree()
            mt = self._build_test_tree(h_lib, twoways_new)
            mt.write_pdbs("alt")
            print len(mt.nodes)


            exit()

            for i, n in enumerate(mt.nodes):
                if i == 0:
                    continue

                print n.motif.name

                if n.motif.mtype == rnamake.motif_type.HELIX:
                    node = mt2.add_motif(h_lib.get_motif(n.motif.name))
                    if node is None:
                        print "fail"
                        exit()

                else:
                    name = n.motif.name
                    try:
                        m1 = twoways_new.get_motif(name+"-0")
                    except:
                        m1 = None

                    try:
                        m2 = twoways_new.get_motif(name+"-1")
                    except:
                        m2 = None

                    node = None

                    fail = 0
                    if m1 is not None:
                        node = mt2.add_motif(m1)
                        if node is not None:
                            end = node.available_ends()[0]
                            ei = node.motif.ends.index(end)

                            d1 = end.d()
                            d2 = n.motif.ends[ei].d()
                            d =rnamake.util.distance(d1, d2)
                            if d > 0.1:
                                fail = 1
                                mt2.remove_node(mt2.last_node)

                        else:
                            fail = 1

                    if not fail:
                        continue

                    fail =0

                    if m2 is not None:
                        node = mt2.add_motif(m2)
                        if node is not None:
                            end = node.available_ends()[0]
                            ei = node.motif.ends.index(end)

                            d1 = end.d()
                            d2 = n.motif.ends[ei].d()
                            d =rnamake.util.distance(d1, d2)
                            if d > 0.1:
                                fail = 1
                                mt2.remove_node(mt2.last_node)
                        else:
                            print "fail to add"
                            exit()



                    if fail:
                        print "stuck"
                        mt.write_pdbs("org")
                        mt2.write_pdbs()
                        exit()


            mt.write_pdbs("org")
            mt2.write_pdbs()
            exit()






def main():
    unittest.main()

if __name__ == '__main__':
    main()
