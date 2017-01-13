import unittest
from rnamake import transform, transformations, sqlite_library
from rnamake.unittests import util
from rnamake.unittests.instances.transform_instances import transform_random
import numpy as np
import random


class BasepairStateAlignmentUnittests(unittest.TestCase):
    """
    Tries to break basepair state alignment. It is critical that this behaves
    properly for MotifStateSearch
    """

    def setUp(self):
        mlib = sqlite_library.MotifSqliteLibrary("ideal_helices")
        self.m = mlib.get(name="HELIX.IDEAL.1")

    def test_align(self):
        m = self.m
        end1 = m.get_end(0).get_state()
        end2 = m.get_end(1).get_state()

        for i in range(100):
            r, t = end1.get_transforming_r_and_t_w_state(end2)
            t += end1.d

            trans = transform.Transform(r, t)
            end2.transform(trans)

            diff = end1.diff(end2)

            if diff > 0.001:
                self.fail("did not align properly")


        end1 = m.get_end(0).get_state()
        end2 = m.get_end(1).get_state()

        for i in range(100):
            r, t = end1.get_transforming_r_and_t_w_state(end2)
            t += end1.d

            trans = transform.Transform(r, t)
            end2.transform(trans)

            diff = end1.diff(end2)
            if diff > 0.01:
                self.fail("did not align properly")

    def test_align_2(self):
        m = self.m
        end1 = m.get_end(0).get_state()
        end2 = m.get_end(1).get_state()

        for i in range(100):
            end2.transform(transform_random())

            r, t = end1.get_transforming_r_and_t_w_state(end2)
            t += end1.d

            trans = transform.Transform(r, t)
            end2.transform(trans)

            diff = end1.diff(end2)
            if diff > 0.001:
                self.fail("did not align properly")

    def test_align_3(self):
        m = self.m
        end1 = m.get_end(0).get_state()
        end2 = m.get_end(0).get_state()

        for i in range(100):
            r, t = end1.get_transforming_r_and_t_w_state(end2)
            t += end1.d

            trans = transform.Transform(r, t)
            end2.transform(trans)

            diff = end1.diff(end2)
            if diff > 0.001:
                self.fail("did not align properly")

    def test_align_4(self):
        m = instances.motif()
        end1 = m.ends[0].state()
        end2 = instances.basepairstate_random()

        for i in range(100):
            r, t = end1.get_transforming_r_and_t_w_state(end2)
            t += end1.d

            new_r, new_d, new_sug = end2.get_transformed_state(r, t)
            end2.set(new_r, new_d, new_sug)

            diff = end1.diff(end2)
            if diff > 0.001:
                self.fail("did not align properly")

    def test_align_5(self):

        for i in range(100):
            m = instances.motif()
            end1 = instances.basepairstate_random()
            end2 = instances.basepairstate_random()

            r, t = end1.get_transforming_r_and_t_w_state(end2)
            t += end1.d

            new_r, new_d, new_sug = end2.get_transformed_state(r, t)
            end2.set(new_r, new_d, new_sug)

            diff = end1.diff(end2)
            if diff > 0.001:
                self.fail("did not align properly")


def main():
    unittest.main()

if __name__ == '__main__':
    main()
