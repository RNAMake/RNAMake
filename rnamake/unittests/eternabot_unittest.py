import unittest
import rnamake.eternabot.sequence_designer as sequence_designer

class EternaBotUnittest(unittest.TestCase):

    def test_creation(self):
        designer = sequence_designer.SequenceDesigner()
        solutions = designer.design("(((....)))", "NNNNNNNNNN")
        print solutions

def main():
    unittest.main()

if __name__ == '__main__':
    main()
