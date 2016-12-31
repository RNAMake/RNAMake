import unittest
import numpy as np

from rnamake import motif_graph, util, exceptions
from rnamake import motif_topology, settings
from rnamake import resource_manager as rm
import build
import secondary_structure_tools
import is_equal


class MotifGraphUnittest(unittest.TestCase):

    def test_creation(self):
        mg = motif_graph.MotifGraph()

    def test_add_motif(self):
        mg = motif_graph.MotifGraph()
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        mg.add_motif(m1)

        # can never use parent_end_index=0 for a graph as that is where that node
        # is already connected to another node
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_motif(m2, parent_end_index=0)

        # supplied parent_end_index and parent_end_name
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_motif(m2, parent_end_index=1, parent_end_name="A1-A8")

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
            mg.add_motif(m2, parent_end_name="A4-A5")

        # invalid parent_end_name, cannot be found as an end in motif
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_motif(m2, parent_end_name="FAKE")

    def test_remove(self):
        builder = build.BuildMotifTree()
        mt = builder.build(2)
        mg = motif_graph.MotifGraph()


        mg2 = motif_graph.MotifGraph()
        mg2.option('sterics', 0)

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        mg.remove_motif(1)

        if len(mg) != 1:
            self.fail("did not remove motif correctly")

    def test_merge(self):
        mg = motif_graph.MotifGraph()
        mg.add_motif(m_name="HELIX.IDEAL")
        mg.add_motif(m_name="HELIX.IDEAL")
        rna_struc = mg.get_structure()

        self.failUnless(len(rna_struc.chains()) == 2)
        self.failUnless(len(rna_struc.residues()) == 6)

    def test_add_connection(self):
        # this is not a valid test as this topology cannot exist!
        mg = motif_graph.MotifGraph()
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m3 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        nway = rm.manager.get_motif(name="NWAY.1GID.0")
        mg.add_motif(m1)
        mg.add_motif(nway)
        mg.add_motif(m2)

        # try connecting through 0th end position
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_connection(1, 2, "A138-A180")

        # try connecting thru an already used end position
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_connection(1, 2, "A141-A162")

        mg.add_connection(1, 2)
        rna_struc = mg.get_structure()
        self.failUnless(len(rna_struc.chains()) == 1)
        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_motif(m3, parent_end_index=1)

        self.failUnless(mg.add_motif(m3) == -1)

        with self.assertRaises(exceptions.MotifGraphException):
            mg.add_connection(1, 2)

    def test_add_connection_2(self):
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        tc = rm.manager.get_motif(name="TC.1S72.0")

        mg = motif_graph.MotifGraph()
        mg.add_motif(tc)
        mg.add_motif(m1)
        mg.add_motif(m2, parent_index=0)
        mg_copy = mg.copy()
        mg.add_connection(0, 2)

        rna_struc_1 = mg.get_structure()
        rna_struc_2 = mg_copy.get_structure()

        self.failUnless(len(rna_struc_1.residues()) == len(rna_struc_2.residues())-2)

    def test_copy(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        new_mg = mg.copy()

        if len(mg.graph) != len(new_mg.graph):
            self.fail("did not copy correctly")

        new_mg.remove_motif(1)

    def test_replace_ideal_helices(self):
        builder = build.BuildMotifGraph()
        mg = builder.build(3)

        new_mg = mg.copy()
        new_mg.replace_ideal_helices()

        struc1 = mg.get_structure()
        struc2 = new_mg.get_structure()

        d1 = struc1.ends[0].d()
        ds_2 = [ struc2.ends[0].d(), struc2.ends[1].d()]

        dist_1 = util.distance(d1, ds_2[0])
        dist_2 = util.distance(d1, ds_2[1])

        self.failIf(dist_1 > 1 and dist_2 > 1,
                    "replacing ideal helices messed up graph")

        if len(new_mg.merger.get_structure().chains()) != 2:
            self.fail("does not have the right number of chains")

    def test_replace_ideal_helices_2(self):
        # specifically test if non aligned nodes are not at the default alignment
        # point that the first HELIX.IDEAL is placed in the right place
        mg = motif_graph.MotifGraph()
        m = rm.manager.get_motif(name='HELIX.IDEAL.6')
        m.move(np.array([40,0,0]))
        mg.add_motif(m)

        new_mg = mg.copy()
        new_mg.replace_ideal_helices()

        struc1 = mg.get_structure()
        struc2 = new_mg.get_structure()

        d1 = struc1.ends[0].d()
        ds_2 = [ struc2.ends[0].d(), struc2.ends[1].d()]

        dist_1 = util.distance(d1, ds_2[0])
        dist_2 = util.distance(d1, ds_2[1])

        self.failIf(dist_1 > 1 and dist_2 > 1,
                    "replacing ideal helices messed up graph")

    def test_replace_helical_sequence(self):
        builder = build.BuildMotifGraph()
        mg = builder.build(3)
        mg.replace_ideal_helices()
        dss = mg.designable_secondary_structure()
        build.fill_basepairs_in_ss(dss)
        mg.replace_helix_sequence(dss)

        self.failUnless(mg.sequence() == dss.sequence())

    def test_replace_helical_sequence_2(self):
        mg = motif_graph.MotifGraph()
        m = rm.manager.get_motif(name='HELIX.IDEAL.6')
        m.move(np.array([40,0,0]))
        mg.add_motif(m)
        mg.replace_ideal_helices()

        new_mg = mg.copy()
        dss = new_mg.designable_secondary_structure()
        build.fill_basepairs_in_ss(dss)
        new_mg.replace_helix_sequence(dss)

        struc1 = mg.get_structure()
        struc2 = new_mg.get_structure()

        d1 = struc1.ends[0].d()
        ds_2 = [ struc2.ends[0].d(), struc2.ends[1].d()]

        dist_1 = util.distance(d1, ds_2[0])
        dist_2 = util.distance(d1, ds_2[1])

        self.failIf(dist_1 > 10 and dist_2 > 10,
                    "replacing helix sequences messed up graph")

    def test_secondary_structure(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()
        for n in mt.tree.nodes:
            mg.add_motif(n.data)
        #mg.write_pdbs("org")
        mg.replace_ideal_helices()

        ss = mg.designable_secondary_structure()

    def test_designable_secondary_structure(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        mg.replace_ideal_helices()
        ss = mg.designable_secondary_structure()

    def test_designable_secondary_structure_2(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        hairpin = rm.manager.mlibs['hairpin'].get_random()
        mg.add_motif(hairpin)

    def test_to_tree(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        mt2 = motif_topology.graph_to_tree(mg)

    def test_mg_to_str(self):
        #load motif graph with all atoms
        base_dir = settings.UNITTEST_PATH + "resources/motif_graph/"
        f = open(base_dir+"test.mg")
        l = f.readline()
        f.close()

        mg = motif_graph.MotifGraph(mg_str=l)

    def test_mg_to_str_multiple_alignments(self):
        base_dir = settings.UNITTEST_PATH + "resources/motif_graph/"
        f = open(base_dir+"tecto_chip_only.mg")
        l = f.readline()
        f.close()

        mg = motif_graph.MotifGraph(mg_str=l)
        self.failIf(len(mg.secondary_structure().chains()) != 1)

    def test_replace_motif(self):
        base_dir = settings.UNITTEST_PATH + "resources/motif_graph/"
        f = open(base_dir+"tecto_chip_only.mg")
        l = f.readline()
        f.close()

        mg = motif_graph.MotifGraph(mg_str=l)
        org_sequence = mg.sequence()
        org_name = mg.get_node(14).data.name
        org_end_name = mg.get_node(14).data.ends[0].name()

        mg.set_options('sterics', 0)
        mg.replace_motif(14, rm.manager.get_motif(name='HELIX.IDEAL.6'))
        self.failUnless(mg.get_node(14).data.name == 'HELIX.IDEAL.6')
        self.failUnless(len(mg.secondary_structure().chains()) == 1)
        self.failIf(mg.sequence() == org_sequence)

        mg.replace_motif(14, rm.manager.get_motif(name=org_name, end_name=org_end_name))
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

    def test_pretty_str(self):
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

    def _test_ss_convert(self):
        mg = motif_graph.MotifGraph()
        mg.add_motif(m_name="HELIX.IDEAL.6")
        mg.add_motif(m_name="TWOWAY.1GID.12")
        mg.add_motif(m_name="HELIX.IDEAL.6")
        mg.replace_ideal_helices()

        for n in mg:
            m = n.data
            name = m.ends[0].res1.name+m.ends[0].res2.name+"="
            name += m.ends[1].res1.name+m.ends[1].res2.name
            print name, m.sequence()

        ss = mg.secondary_structure()
        print "made it"
        for m in ss.motifs:
            name = m.ends[0].res1.name+m.ends[0].res2.name+"="
            name += m.ends[1].res1.name+m.ends[1].res2.name
            print name, m.sequence()


def main():
    unittest.main()

if __name__ == '__main__':
    main()














