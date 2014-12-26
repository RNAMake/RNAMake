import unittest
import rnamake.motif_tree_state

class MotifTreeStateUnittest(unittest.TestCase):

    def test_creation(self):
        pass

    def test_parse_db_name(self):
        name = "HELIX.LE.16-0-0-0-0-1-1"
        name_elements = rnamake.motif_tree_state.parse_db_name(name)
        if name_elements.motif_name != "HELIX.LE.16":
            self.fail()

def main():
    unittest.main()

if __name__ == '__main__':
    main()
