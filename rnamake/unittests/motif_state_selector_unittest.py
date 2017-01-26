import unittest

from rnamake import resource_manager, motif_state_selector, motif_state


class MotifStateSelectorUnittest(unittest.TestCase):

    def setUp(self):
        self.rm = resource_manager.ResourceManager()

    def test_add(self):
        selector = motif_state_selector.MotifStateSelector(self.rm)
        selector.add('ideal_helices')
        connected = selector.get_connected_libs(-1)
        self.failUnless(len(connected) == 1)
        self.failUnless(connected[0] == 0)
        motif_states = selector.get_motif_states(connected[0])
        self.failUnless(len(motif_states) > 0)

    def test_add_2(self):
        selector = motif_state_selector.MotifStateSelector(self.rm)
        selector.add('unique_twoway')
        selector.add('ideal_helices_min')

        with self.assertRaises(ValueError):
            selector.add('FAKE')

    def test_add_3(self):
        m1 = self.rm.get_state(name='HELIX.IDEAL')
        m2 = self.rm.get_state(name='TWOWAY.2VQE.19', end_name='A1008-A1021')

        mlib1 = motif_state.MotifLibrary('helix', [m1])
        mlib2 = motif_state.MotifLibrary('twoway', [m2])

        selector = motif_state_selector.MotifStateSelector(self.rm)
        selector.add(mlib=mlib1)
        selector.add(mlib=mlib2)
        selector.connect('helix', 'twoway')

        self.failUnless(selector.get_connected_libs(0)[0] == 1)
        self.failUnless(selector.get_connected_libs(1)[0] == 0)

    def test_default_selector(self):
        selector = motif_state_selector.default_selector(self.rm)



def main():
    unittest.main()

if __name__ == '__main__':
    main()