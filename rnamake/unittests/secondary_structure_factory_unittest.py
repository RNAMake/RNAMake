import unittest
import rnamake.secondary_structure_factory as secondary_structure_factory
import rnamake.resource_manager as rm

class SecondaryStructureFactoryUnittest(unittest.TestCase):

    def test_creation(self):
        m = rm.manager.get_motif("NWAY.1S72.18-01551-01634")
        seq = "GCAU+AUGC"
        db  = "((((+))))"

        print seq, db
        ss = secondary_structure_factory.factory.get_structure(seq, db)



def main():
    unittest.main()

if __name__ == '__main__':
    main()