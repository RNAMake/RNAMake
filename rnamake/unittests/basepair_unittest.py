import unittest
from rnamake import structure, transformations, x3dna, residue_type, basepair, settings

import util, instances
import numpy as np
import copy
import random

from instances.transform_instances import transform_random
import numerical

def _calc_center(res):
    center = np.array([0.0,0.0,0.0])
    count = 0
    for r in res:
        for a in r:
            if a is None:
                continue
            center += a.coords
            count += 1
    center /= count
    return center

def _calc_name(res):
        res1, res2 = res

        res1_name = res1.chain_id+str(res1.num)+str(res1.i_code)
        res2_name = res2.chain_id+str(res2.num)+str(res2.i_code)

        if res1.chain_id < res2.chain_id:
            return res1_name+"-"+res2_name
        if res1.chain_id > res2.chain_id:
            return res2_name+"-"+res1_name

        if res1.num < res2.num:
            return res1_name+"-"+res2_name
        else:
            return res2_name+"-"+res1_name


class BasepairUnittest(unittest.TestCase):

    def setUp(self):
        self.rts = residue_type.ResidueTypeSet()
        path = settings.UNITTEST_PATH + "resources/motifs/p4p6/p4p6.pdb"
        s = structure.structure_from_pdb(path, self.rts)

        x = x3dna.X3dna()
        x_bps = x.get_basepairs(settings.UNITTEST_PATH + "resources/motifs/p4p6/p4p6.pdb")
        x_bp = x_bps[0]
        res1 = s.get_residue(num=x_bp.res1.num)
        res2 = s.get_residue(num=x_bp.res2.num)
        res = [res1, res2]
        center = _calc_center(res)
        bp = basepair.Basepair(res1.uuid, res2.uuid, x_bp.r, center,
                               [res1.get_coords("C1'"), res2.get_coords("C1'")],
                               _calc_name(res), bp_type=x_bp.bp_type)

        self.basepair = bp

    def test_synced(self):
        rts = residue_type.ResidueTypeSet()
        path = settings.UNITTEST_PATH + "resources/motifs/p4p6/p4p6.pdb"
        s = structure.structure_from_pdb(path, rts)

        x = x3dna.X3dna()
        x_bps = x.get_basepairs(settings.UNITTEST_PATH + "resources/motifs/p4p6/p4p6.pdb")
        x_bp = x_bps[0]
        res1 = s.get_residue(num=x_bp.res1.num)
        res2 = s.get_residue(num=x_bp.res2.num)
        res = [res1, res2]
        center = _calc_center(res)
        bp = basepair.Basepair(res1.uuid, res2.uuid, x_bp.r, center,
                               [res1.get_coords("C1'"), res2.get_coords("C1'")],
                               _calc_name(res), bp_type=x_bp.bp_type)


        for i in range(100):
            t = transform_random()

            bp.transform(t)

            for r in [res1, res2]:
                r.transform(t)

        self.failUnless(numerical.are_points_equal(bp.d, _calc_center(res)))
        self.failUnless(numerical.are_points_equal(bp.res1_sugar, res1.get_coords("C1'")))
        self.failUnless(numerical.are_points_equal(bp.res2_sugar, res2.get_coords("C1'")))

    def test_copy(self):
        bp = self.basepair
        bp_copy = basepair.Basepair.copy(bp)


"""class BasepairStateUnittest(unittest.TestCase):

    def test_copy(self):
        bpstate = rnamake.basepair.BasepairState(np.eye(3),
                                                 np.array([1, 0, 0]),
                                                 [[1, 0, 0], [0, 1, 0]])

        cbpstate = bpstate.copy()
        cbpstate.r[0][0] += 1
        if bpstate.r[0][0] == cbpstate.r[0][0]:
            self.fail()

        cbpstate.sugars[0][0] += 1
        if bpstate.sugars[0][0] == cbpstate.sugars[0][0]:
            self.fail("sugars did not deep copy")

"""

def main():
    unittest.main()

if __name__ == '__main__':
    main()
