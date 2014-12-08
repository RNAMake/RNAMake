import unittest
import rnamake.pdb_parser
import rnamake

class PdbParserUnittest(unittest.TestCase):

    def test_parse(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6.pdb"
        residues = rnamake.pdb_parser.parse(path)

def main():
    unittest.main()

if __name__ == '__main__':
    main()
