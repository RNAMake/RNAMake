import unittest
import os
from rnamake import settings
from rnamake.wrappers import simulate_tectos_wrapper
from rnamake import motif_graph



RESOURCE_PATH = settings.LIB_PATH + "/lib/RNAMake/unittests/unittest_resources/apps/design_rna/"


class SimulateTectosUnittest(unittest.TestCase):

    def test_wildtype(self):
        stw = simulate_tectos_wrapper.SimulateTectosWrapper()
        stw.run()
        hit_count = stw.get_output()
        self.failUnless(hit_count > 970 and hit_count < 1250,
                        "wildtype did not have expected hit_count: " + str(hit_count))

    def test_high(self):
        stw = simulate_tectos_wrapper.SimulateTectosWrapper()
        stw.run(cseq="CTAGGATATGGCTCTAGGGGGGGAACCCCCTAGAGCCTAAGTCCTAG")
        hit_count = stw.get_output()
        self.failUnless(hit_count > 6900 and hit_count < 7750,
                        "best seq did not have expected hit_count: " + str(hit_count))

    def test_low(self):
        stw = simulate_tectos_wrapper.SimulateTectosWrapper()
        stw.run(cseq="CTAGGATATGGAGACGCCTGGGGAACCAGGCGTCTCCTAAGTCCTAG")
        hit_count = stw.get_output()
        self.failUnless(hit_count > 800 and hit_count < 1250,
                        "worst seq did not have expected hit_count: " + str(hit_count))

def main():
    unittest.main()

if __name__ == '__main__':
    main()

