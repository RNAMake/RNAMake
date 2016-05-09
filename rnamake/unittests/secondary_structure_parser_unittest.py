import unittest
import glob

from rnamake import secondary_structure_parser, settings, util, exceptions

def add_new_testcase(seq, db, g):
    path = settings.UNITTEST_PATH + "/resources/secondary_structure/"
    files = glob.glob(path+"/*")
    pos =  len(files)
    f = open(path + "/case." + str(pos) + ".dat", "w")
    f.write(seq + "\n" + db + "\n" + str(g))
    f.close()


class SecondaryStructureChainGraph(unittest.TestCase):
    def test_creation(self):
        secondary_structure_parser.SecondaryStructureChainGraph()


class SecondaryStructureParserUnittest(unittest.TestCase):

    def test_creation(self):
        p = secondary_structure_parser.SecondaryStructureParser()

        with self.assertRaises(exceptions.SecondaryStructureParserException):
            p.parse()

        with self.assertRaises(exceptions.SecondaryStructureParserException):
            p.parse("GG+CC", "(((+))")

        with self.assertRaises(exceptions.SecondaryStructureParserException):
            p.parse("GG+CC")

        with self.assertRaises(exceptions.SecondaryStructureParserException):
            p.parse(dot_bracket="((+))")

    def test_parse(self):
        p = secondary_structure_parser.SecondaryStructureParser()

        with self.assertRaises(exceptions.SecondaryStructureParserException):
            p.parse("GG+CC", "()+))")

        p.reset()
        g = p.parse("GGG+CC", ".((+))")
        self.failUnless(len(g) == 5, "did not parse left unpaired correctly")

        p.reset()
        g = p.parse("GG+CCAA", "((+))..")
        # 5 nodes not 5 residues
        self.failUnless(len(g) == 5, "did not parse left unpaired correctly")


    def test_parse_to_motifs(self):
        parser = secondary_structure_parser.SecondaryStructureParser()
        seq = "GGAAGACAAGACAACC"
        ss  = "((..(.)..(.)..))"
        #seq = "GG+CC"
        #ss  = "((+))"
        motifs = parser.parse_to_motifs(seq, ss)
        #print motifs
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