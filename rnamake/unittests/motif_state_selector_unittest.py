import unittest
import build
import rnamake.motif_state_search
import rnamake.motif_state_selector

class MotifStateSelectorUnittest(unittest.TestCase):

    def _get_dummy_node(self):
        builder = build.BuildMotifTree()
        mt = builder.build()
        start = mt.get_node(0).data.ends[0].state()
        end   = mt.get_node(1).data.ends[1].state()
        mss = rnamake.motif_state_search.MotifStateSearch()
        n = mss._start_node(start)
        return n

    def test_creation(self):
        selector = rnamake.motif_state_selector.MotifStateSelector()
        selector.add('twoway')
        selector.connect('twoway', 'twoway')

    def test_get_children_ms(self):
        n = self._get_dummy_node()
        selector = rnamake.motif_state_selector.default_selector()
        motif_states, types = selector.get_children_ms(n)
        if len(motif_states) == 0:
            self.fail("did not get any motif_states back")

        n.ntype = 0
        motif_states, types = selector.get_children_ms(n)
        if len(motif_states) == 0:
            self.fail("did not get any motif_states back")

    def test_round_robin(self):
        n = self._get_dummy_node()
        selector = rnamake.motif_state_selector.MSS_RoundRobin()
        selector.add("twoway")
        selector.add("ideal_helices")

        n.ntype = 0
        motif_states, types = selector.get_children_ms(n)
        unique_types = []
        for t in types:
            if t not in unique_types:
                unique_types.append(t)
        if len(unique_types) != 2:
            self.fail("did not get the correct motif_states back")

    def test_helix_flank(self):
        n = self._get_dummy_node()
        selector = rnamake.motif_state_selector.MSS_HelixFlank()
        selector.add("twoway")

        motif_states, types = selector.get_children_ms(n)

        n.ntype = 0
        motif_states, types = selector.get_children_ms(n)
        if len(motif_states) < 100:
            self.fail("did not get the expected number of states")

def main():
    unittest.main()

if __name__ == '__main__':
    main()