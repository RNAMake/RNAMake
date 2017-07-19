import unittest

from rnamake import motif, sqlite_library

class MotifAlignerUnittest(unittest.TestCase):

    def test_creation(self):
        pass

    def test_align(self):
        mlib = sqlite_library.MotifSqliteLibrary('ideal_helices')
        m1 = mlib.get(name='HELIX.IDEAL.2')
        m2 = mlib.get(name='HELIX.IDEAL.2')
        m3 = mlib.get(name='HELIX.IDEAL.2')

        m_aligner = motif.MotifAligner()
        m_aligner.align(m1.get_end(1), m2)
        m_aligner.align(m2.get_end(1), m3)

        diff = m1.get_end(1).diff(m2.get_end(0))
        self.failUnless(diff < 0.001)

        diff = m2.get_end(1).diff(m3.get_end(0))
        self.failUnless(diff < 0.001)

    def test_get_aligned(self):
        mlib = sqlite_library.MotifSqliteLibrary('ideal_helices')
        m1 = mlib.get(name='HELIX.IDEAL.2')
        m2 = mlib.get(name='HELIX.IDEAL.2')

        m_aligner = motif.MotifAligner()
        m_aligned = m_aligner.get_aligned(m1.get_end(1), m2)

        diff = m1.get_end(1).diff(m_aligned.get_end(0))
        self.failUnless(diff < 0.001)

        diff = m1.get_end(1).diff(m2.get_end(0))
        self.failUnless(diff > 0.001)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
