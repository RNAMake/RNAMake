import unittest
from rnamake import transformations
from rnamake.unittests import util, instances
import numpy as np
import random


class BasepairStateAlignmentUnittests(unittest.TestCase):
    """
    Tries to break basepair state alignment. It is critical that this behaves
    properly for MotifStateSearch
    """

    def test_align(self):
        m = instances.motif()
        end1 = m.ends[0].state().copy()
        end2 = m.ends[1].state().copy()


        for i in range(100):
            r, t = end1.get_transforming_r_and_t_w_state(end2)
            t += end1.d

            new_r, new_d, new_sug = end2.get_transformed_state(r, t)
            end2.set(new_r, new_d, new_sug)

            diff = end1.diff(end2)

            print end1.sugars
            print end2.sugars
            exit()
            if diff > 0.001:
                self.fail("did not align properly")

        end1 = m.ends[1].state().copy()
        end2 = m.ends[0].state().copy()

        for i in range(100):
            r, t = end1.get_transforming_r_and_t_w_state(end2)
            t += end1.d

            new_r, new_d, new_sug = end2.get_transformed_state(r, t)
            end2.set(new_r, new_d, new_sug)

            diff = end1.diff(end2)
            if diff > 0.01:
                self.fail("did not align properly")

    def test_align_2(self):
        m = instances.motif()
        end1 = m.ends[0].state()
        end2 = m.ends[1].state()

        for i in range(100):
            trans = transformations.random_rotation_matrix()
            rand_r = trans[:3,:3]
            rand_t = np.array([random.randint(0, 10),
                               random.randint(0, 10),
                               random.randint(0, 10)])
            new_r, new_d, new_sug = end2.get_transformed_state(rand_r, rand_t)
            end2.set(new_r, new_d, new_sug)

            r, t = end1.get_transforming_r_and_t_w_state(end2)
            t += end1.d

            new_r, new_d, new_sug = end2.get_transformed_state(r, t)
            end2.set(new_r, new_d, new_sug)

            diff = end1.diff(end2)
            if diff > 0.001:
                self.fail("did not align properly")

    def test_align_3(self):
        m = instances.motif()
        end1 = m.ends[0].state()
        end2 = m.ends[0].state()

        for i in range(100):
            r, t = end1.get_transforming_r_and_t_w_state(end2)
            t += end1.d

            new_r, new_d, new_sug = end2.get_transformed_state(r, t)
            end2.set(new_r, new_d, new_sug)

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
