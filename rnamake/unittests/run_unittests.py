import unittest

testmodules = [
    'atom_unittest',
    'base_unittest',
    'basepair_unittest',
    'basic_io_unittest',
    'chain_unittest',
    'eternabot_unittest',
    'graph_unittest',
    'motif_unittest',
    'motif_graph_unittest',
    'motif_library_unittest',
    'motif_scorer_unittest',
    'motif_state_search_unittest',
    'motif_state_selector_unittest',
    'motif_state_ensemble_tree_unittest',
    'motif_state_tree_unittest',
    'motif_tree_unittest',
    'option_unittest',
    'pdb_parser_unittest',
    'pose_unittest',
    'residue_unittest',
    'residue_type_unittest',
    'resource_manager_unittest',
    'secondary_structure_unittest',
    'segmenter_unittest',
    'sqlite_library_unittest',
    'steric_lookup_unittest',
    'structure_unittest',
    'transform_unittest',
    'thermo_fluc_sampler_unittest',
    'tree_unittest',
    'vienna_unittest',
    'x3dna_unittest'
    ]

suite = unittest.TestSuite()

for t in testmodules:
    try:
        # If the module defines a suite() function, call it to get the suite.
        mod = __import__(t, globals(), locals(), ['suite'])
        suitefn = getattr(mod, 'suite')
        suite.addTest(suitefn())
    except (ImportError, AttributeError):
        # else, just load all the test cases from the module.
        suite.addTest(unittest.defaultTestLoader.loadTestsFromName(t))

unittest.TextTestRunner().run(suite)
