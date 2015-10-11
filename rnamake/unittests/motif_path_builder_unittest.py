import unittest
import rnamake.motif_path_builder as motif_path_builder

class MotifPathBuilderUnittest(unittest.TestCase):

    def test_creation(self):
        builder = motif_path_builder.MotifPathBuilder()

def main():
    unittest.main()

if __name__ == '__main__':
    main()