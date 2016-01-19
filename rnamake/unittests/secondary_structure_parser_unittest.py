import unittest
import glob
from rnamake import secondary_structure_parser, settings, util

def add_new_testcase(seq, db, g):
    path = settings.UNITTEST_PATH + "/resources/secondary_structure/"
    files = glob.glob(path+"/*")
    pos =  len(files)
    f = open(path + "/case." + str(pos) + ".dat", "w")
    f.write(seq + "\n" + db + "\n" + str(g))
    f.close()


class SecondaryStructureParserUnittest(unittest.TestCase):

    def test_creation(self):
        parser = secondary_structure_parser.SecondaryStructureParser()

    def _test_parse(self):
        path = settings.UNITTEST_PATH + "/resources/secondary_structure/"
        files = glob.glob(path+"/*")
        for f_name in files:
            f = open(f_name)
            lines = f.readlines()
            f.close()

            name = util.filename(f_name)

            seq, db, expected_g = lines[0].rstrip(), lines[1].rstrip(), lines[2:]
            expected_g = "".join(expected_g)
            parser = secondary_structure_parser.SecondaryStructureParser()
            g = parser.parse(seq, db)
            if str(g) != expected_g:
                self.fail("failed case: " + name + " in test parse")

    def test_parse_new(self):
        parser = secondary_structure_parser.SecondaryStructureParser()
        #seq = "GGAAGACAAGACAACC"
        #ss  = "((..(.)..(.)..))"
        seq = "GG+CC"
        ss  = "((+))"

        g = parser.parse(seq, ss)
        print g

        #add_new_testcase(seq, ss, g)

    def test_parse_to_motifs(self):
        parser = secondary_structure_parser.SecondaryStructureParser()
        #seq = "GGAAGACAAGACAACC"
        #ss  = "((..(.)..(.)..))"
        seq = "GG+CC"
        ss  = "((+))"
        motifs = parser.parse_to_motifs(seq, ss)
        print motifs
        #if len(motifs) != 4:
        #    self.fail("did not get the right number of motifs")

    def test_parse_to_motif_graph(self):
        parser = secondary_structure_parser.SecondaryStructureParser()
        seq = "GGAAGACAAGACAACC"
        ss  = "((..(.)..(.)..))"
        g = parser.parse_to_motif_graph(seq, ss)
        if len(g.graph) != 4:
            self.fail("did not build graph properly")


def main():
    unittest.main()

if __name__ == '__main__':
    main()