import unittest
import simulate_tectos

from rnamake import resource_manager as rm
from rnamake import sqlite_library, basepair

class ThermoSamplerTests(unittest.TestCase):

    def test_simulate_tectos_start(self):
        st = simulate_tectos.SimulateTectos()
        mst = st.mset.to_mst()
        first_end = mst.get_node(1).data.cur_state.end_states[1]
        last_end = mst.last_node().data.cur_state.end_states[1]

        diff = first_end.diff(last_end)
        expected_diff = 10.9928543144
        self.failIf(abs(diff - expected_diff) > 0.1)




def main():
    unittest.main()

if __name__ == '__main__':
    main()
