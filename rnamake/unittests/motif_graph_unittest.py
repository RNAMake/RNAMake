import unittest
import numpy as np

from rnamake import motif_graph, util, exceptions
from rnamake import motif_topology, settings
from rnamake import resource_manager, sqlite_library
import build
import secondary_structure_tools
import is_equal


class MotifGraphUnittest(unittest.TestCase):

    def setUp(self):
        self.rm = resource_manager.ResourceManager()

    def test_add_motif(self):
        mg = motif_graph.MotifGraph(self.rm)
        m1 = self.rm.get_motif(name="HELIX.IDEAL.2")
        m2 = self.rm.get_motif(name="HELIX.IDEAL.2")
        mg.add_motif(m1)

        # can never use parent_end_index=0 for a graph as that is where that node
        # is already connected to another node
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_motif(m2, parent_end_index=0)

        # supplied parent_end_index and parent_end_name
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_motif(m2, parent_end_index=1, parent_end_name="A4-A5")

        # must supply a motif or motif name
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_motif()

        # motif not found in resource manager
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_motif(m_name="FAKE")

        # catches invalid parent_index
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_motif(m2, parent_index=2)

        # invalid parent_end_index, has only 0 and 1
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_motif(m2, parent_end_index=3)

        # invalid parent_end_name, is the name of end 0
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_motif(m2, parent_end_name="A1-A8")

        # invalid parent_end_name, cannot be found as an end in motif
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_motif(m2, parent_end_name="FAKE")

    def test_remove(self):
        builder = build.BuildMotifTree(self.rm)
        mt = builder.build(2)
        mg = motif_graph.MotifGraph(self.rm)

        for n in mt:
            mg.add_motif(n.data)

        mg.remove_motif(1)

        if len(mg) != 1:
            self.fail("did not remove motif correctly")

    def test_merge(self):
        mg = motif_graph.MotifGraph(self.rm)
        mg.add_motif(m_name="HELIX.IDEAL")
        mg.add_motif(m_name="HELIX.IDEAL")
        rna_struc = mg.get_structure()

        self.failUnless(rna_struc.num_chains() == 2)
        self.failUnless(rna_struc.num_res() == 6)

    def test_add_connection(self):
        # this is not a valid test as this topology cannot exist!
        mg = motif_graph.MotifGraph(self.rm)
        mg.add_motif(m_name="HELIX.IDEAL.2")
        mg.add_motif(m_name="NWAY.1GID.0")
        mg.add_motif(m_name="HELIX.IDEAL.2")
        m3 = self.rm.get_motif(name="HELIX.IDEAL.2")

        # try connecting through 0th end position
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_connection(1, 2, "A138-A180")

        # try connecting thru an already used end position
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_connection(1, 2, "A141-A162")

        mg.add_connection(1, 2)
        rna_struc = mg.get_structure()
        self.failUnless(rna_struc.num_chains() == 1)
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_motif(m3, parent_end_index=1)

        self.failUnless(mg.add_motif(m3) == -1)

        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_connection(1, 2)

    def test_add_connection_2(self):
        m1 = self.rm.get_motif(name="HELIX.IDEAL.2")
        m2 = self.rm.get_motif(name="HELIX.IDEAL.2")
        tc = self.rm.get_motif(name="TC.1S72.0")

        mg = motif_graph.MotifGraph(self.rm)
        mg.add_motif(tc)
        mg.add_motif(m1)
        mg.add_motif(m2, parent_index=0)
        mg_copy =  motif_graph.MotifGraph.copy(mg)
        mg.add_connection(0, 2)
        mg.nodes_to_pdbs()

        rna_struc_1 = mg.get_structure()
        rna_struc_2 = mg_copy.get_structure()

        self.failUnless(rna_struc_1.num_res() == rna_struc_2.num_res()-2)

    def test_copy(self):
        builder = build.BuildMotifTree(self.rm)
        mt = builder.build(3)
        mg = motif_graph.MotifGraph(self.rm)

        for n in mt:
            mg.add_motif(n.data)

        new_mg =  motif_graph.MotifGraph.copy(mg)

        if len(mg) != len(new_mg):
            self.fail("did not copy correctly")

        new_mg.remove_motif(1)

    def test_replace_ideal_helices(self):
        builder = build.BuildMotifGraph(self.rm)
        mg = builder.build(3)

        new_mg = motif_graph.MotifGraph.copy(mg)
        new_mg.replace_ideal_helices()

        struc1 = mg.get_structure()
        struc2 = new_mg.get_structure()

        d1 = struc1.get_end(0).d
        ds_2 = [ struc2.get_end(0).d, struc2.get_end(1).d]

        dist_1 = util.distance(d1, ds_2[0])
        dist_2 = util.distance(d1, ds_2[1])

        self.failIf(dist_1 > 1 and dist_2 > 1,
                    "replacing ideal helices messed up graph")

        rna_struc = new_mg.get_structure()
        if rna_struc.num_chains() != 2:
            self.fail("does not have the right number of chains")

    def test_replace_ideal_helices_2(self):
        # specifically test if non aligned nodes are not at the default alignment
        # point that the first HELIX.IDEAL is placed in the right place
        mg = motif_graph.MotifGraph(self.rm)
        m = self.rm.get_motif(name='HELIX.IDEAL.6')
        m.move(np.array([40,0,0]))
        mg.add_motif(m)

        new_mg = motif_graph.MotifGraph.copy(mg)
        new_mg.replace_ideal_helices()

        struc1 = mg.get_structure()
        struc2 = new_mg.get_structure()

        d1 = struc1.get_end(0).d
        ds_2 = [ struc2.get_end(0).d, struc2.get_end(1).d]

        dist_1 = util.distance(d1, ds_2[0])
        dist_2 = util.distance(d1, ds_2[1])

        self.failIf(dist_1 > 1 and dist_2 > 1,
                    "replacing ideal helices messed up graph")

    def test_bp_steps(self):
        mg = motif_graph.MotifGraph(self.rm)
        mlib = sqlite_library.MotifSqliteLibrary("bp_steps")
        mlib.load_all()
        for m in mlib.all():
            mg.add_motif(m)
        mg.nodes_to_pdbs()

    def test_replace_helical_sequence(self):
        builder = build.BuildMotifGraph(self.rm)
        mg = builder.build(3)
        #for n in mg:
        #    print n.data.name, n.data.get_end(0).name
        mg.replace_ideal_helices()
        dss = mg.designable_secondary_structure()
        build.fill_basepairs_in_ss(dss)
        mg.replace_helix_sequence(dss)

        #print mg.sequence()
        #print dss.sequence()
        #mg.to_pdb("test.pdb", renumber=1, close_chain=1)
        #mg.nodes_to_pdbs()
        self.failUnless(mg.sequence() == dss.sequence())

    def _test_replace_helical_sequence_spec(self):
        self.rm.motif_factory._ref_motif.to_pdb("ref.pdb")
        mg = motif_graph.MotifGraph(self.rm)
        mg.add_motif(m_name="HELIX.IDEAL.15")
        mg.add_motif(m_name="TWOWAY.3R1C.32", m_end_name="g4-h5")
        mg.add_motif(m_name="HELIX.IDEAL.14")
        mg.replace_ideal_helices()
        dss = mg.secondary_structure()
        rna_struc = mg.get_structure()
        mg.to_pdb("org.pdb", renumber=1, close_chain=1)

        #print rna_struc.get_chain(0).get_residue(0).name
        #print dss.get_chain(0).get_residue(0).name
        seq = "ACUGAGGAACGUACGACGGCUUACCGUUUAACUA&UAGUUAAACGGUAAGCGGUCGUACGUUCCUCAGU"
        dss.replace_sequence(seq)
        #print dss.get_chain(0).get_residue(0).name
        mg.replace_helix_sequence(dss)
        print mg.sequence()
        print seq
        rna_struc = mg.get_structure()
        #print rna_struc.get_chain(0).get_residue(0).name
        mg.to_pdb("test.pdb", renumber=1, close_chain=1)
        exit()

        mg.to_pdb("test.pdb", renumber=1, close_chain=1)

    def test_replace_helical_sequence_2(self):
        mg = motif_graph.MotifGraph(self.rm)
        m = self.rm.get_motif(name='HELIX.IDEAL.6')
        m.move(np.array([40,0,0]))
        mg.add_motif(m)
        mg.replace_ideal_helices()

        new_mg = motif_graph.MotifGraph.copy(mg)
        dss = new_mg.designable_secondary_structure()
        build.fill_basepairs_in_ss(dss)

        new_mg.replace_helix_sequence(dss)

        struc1 = mg.get_structure()
        struc2 = new_mg.get_structure()

        d1 = struc1.get_end(0).d
        ds_2 = [ struc2.get_end(0).d, struc2.get_end(1).d]

        dist_1 = util.distance(d1, ds_2[0])
        dist_2 = util.distance(d1, ds_2[1])

        self.failIf(dist_1 > 10 and dist_2 > 10,
                    "replacing helix sequences messed up graph")

    def _test_to_tree(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        mt2 = motif_topology.graph_to_tree(mg)

    def test_mg_to_str(self):
        builder = build.BuildMotifGraph(self.rm)
        mg = builder.build(3)
        s = mg.to_str()

        mg2 = motif_graph.MotifGraph(self.rm, s)
        self.failUnless(len(mg) == len(mg2))

    def _test_mg_to_str(self):
        #load motif graph with all atoms
        base_dir = settings.UNITTEST_PATH + "resources/motif_graph/"
        f = open(base_dir+"test.mg")
        l = f.readline()
        f.close()

        mg = motif_graph.MotifGraph(self.rm, mg_str=l)

    def _test_mg_to_str_multiple_alignments(self):
        base_dir = settings.UNITTEST_PATH + "resources/motif_graph/"
        f = open(base_dir+"tecto_chip_only.mg")
        l = f.readline()
        f.close()

        mg = motif_graph.MotifGraph(mg_str=l)
        self.failIf(len(mg.secondary_structure().chains()) != 1)

    def _test_replace_motif(self):
        base_dir = settings.UNITTEST_PATH + "resources/motif_graph/"
        f = open(base_dir+"tecto_chip_only.mg")
        l = f.readline()
        f.close()

        mg = motif_graph.MotifGraph(self.rm, mg_str=l)
        org_sequence = mg.sequence()
        org_name = mg.get_node(14).data.name
        org_end_name = mg.get_node(14).data.ends[0].name()

        mg.set_sterics(0)
        mg.replace_motif(14, self.rm.get_motif(name='HELIX.IDEAL.6'))
        self.failUnless(mg.get_node(14).data.name == 'HELIX.IDEAL.6')
        self.failUnless(len(mg.secondary_structure().chains()) == 1)
        self.failIf(mg.sequence() == org_sequence)

        mg.replace_motif(14, self.rm.get_motif(name=org_name, end_name=org_end_name))
        self.failUnless(mg.get_node(14).data.name == org_name)
        self.failUnless(mg.sequence() == org_sequence)
        self.failUnless(len(mg.secondary_structure().chains()) == 1)

    def _test_to_str(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        s = mg.to_str()
        mg_new = motif_graph.MotifGraph(mg_str=s)
        if len(mg_new) != len(mg):
            self.fail("did not get right number of nodes")

        atoms1 = mg.get_structure().structure.atoms()
        atoms2 = mg_new.get_structure().structure.atoms()

        for i in range(len(atoms1)):
            if not is_equal.are_atom_equal(atoms1[i], atoms2[i]):
                self.fail("atoms are not equal")

    def _test_pretty_str(self):
        builder = build.BuildMotifGraph()
        mg = builder.build(3)
        mg.to_pretty_str()

    def _test_pretty_str_2(self):
        mg_file = settings.UNITTEST_PATH + "resources/motif_graph/mini_ttr.mg"
        f = open(mg_file)
        lines = f.readlines()
        f.close()
        mg = motif_graph.MotifGraph(mg_str=lines[0])
        #print mg.to_pretty_str()
        #mg.write_pdbs()


def main():
    unittest.main()

if __name__ == '__main__':
    main()














