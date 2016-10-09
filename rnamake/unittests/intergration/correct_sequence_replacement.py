import unittest
from rnamake import sqlite_library, motif_graph, util
from rnamake.unittests import build


class SequenceReplacementUnittest(unittest.TestCase):

    def setUp(self):
        pass

    def _compare_res_pos(self, bp1, bp2):
        dist1 = util.distance(bp1.res1.get_atom("C1'").coords, bp2.res1.get_atom("C1'").coords)
        dist2 = util.distance(bp1.res1.get_atom("C1'").coords, bp2.res2.get_atom("C1'").coords)

        if dist1 < dist2:
            print "made it 1"
            if bp1.res1.name != bp2.res1.name or bp1.res2.name != bp2.res2.name:
                self.fail("do not have the same residue in place")
        else:
            print bp1.res1.name+bp1.res2.name, bp2.res1.name+bp2.res2.name
            if bp1.res1.name != bp2.res2.name or bp1.res2.name != bp2.res1.name:
                self.fail("do not have the same residue in place")


    def test_twoways(self):
        mlib = sqlite_library.MotifSqliteLibrary("twoway")
        mlib.load_all()

        bad_keys = "TWOWAY.2GDI.4-X20-X45 TWOWAY.1S72.46-02097-02647 TWOWAY.2GDI.6-Y20-Y45".split()

        for m in mlib.all():
            #if m.name != "TWOWAY.1S72.29":
            #    continue
            key = m.name+"-"+m.ends[0].name()
            if key in bad_keys:
                continue
            mg = motif_graph.MotifGraph()
            mg.add_motif(m_name="HELIX.IDEAL")
            mg.add_motif(m_name="HELIX.IDEAL")
            mg.add_motif(m_name="HELIX.IDEAL")
            mg.add_motif(m)
            mg.add_motif(m_name="HELIX.IDEAL")
            mg.add_motif(m_name="HELIX.IDEAL")
            mg.add_motif(m_name="HELIX.IDEAL")
            mg.replace_ideal_helices()
            dss = mg.designable_secondary_structure()
            build.fill_basepairs_in_ss(dss)
            mg.replace_helix_sequence(dss)

            print m.name, m.ends[0].name()
            mg.write_pdbs()

            bp1 = mg.get_node(2).data.ends[1]
            bp2 = mg.get_node(3).data.ends[0]

            self._compare_res_pos(bp1, bp2)

            bp1 = mg.get_node(3).data.ends[1]
            bp2 = mg.get_node(4).data.ends[0]

            self._compare_res_pos(bp1, bp2)
            mg.remove_node_level()



def main():
    unittest.main()

if __name__ == '__main__':
    main()