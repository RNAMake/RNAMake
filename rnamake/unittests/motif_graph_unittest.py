import unittest
import rnamake.motif_graph as motif_graph
import rnamake.motif_type as motif_type
import rnamake.graph as graph
import rnamake.util as util
import rnamake.eternabot.sequence_designer as sd
import rnamake.resource_manager as rm
from rnamake import motif_topology, secondary_structure_graph, settings
import build
import numerical
import secondary_structure_tools
import is_equal


class MotifGraphUnittest(unittest.TestCase):

    def test_creation(self):
        mg = motif_graph.MotifGraph()

    def test_add(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)

        mg = motif_graph.MotifGraph()
        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        if len(mg.graph) != 3:
            self.fail("did not get the right number of motifs")

    def test_remove(self):
        builder = build.BuildMotifTree()
        mt = builder.build(2)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        mg.remove_motif(1)

        if len(mg) != 1:
            self.fail("did not remove motif correctly")

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
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        new_mg = mg.copy()
        new_mg.replace_ideal_helices()

        d1 = mg.last_node().data.ends[1].d()
        d2 = new_mg.last_node().data.ends[1].d()
        diff = util.distance(d1, d2)
        if diff > 1:
            self.fail("replacing ideal helices messed up graph")

        if len(new_mg.merger.get_structure().chains()) != 2:
            self.fail("does not have the right number of chains")

    def _test_get_ss(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        new_mg = mg.copy()
        new_mg.replace_ideal_helices()
        ss = new_mg.designable_secondary_structure()
        designer = sd.SequenceDesigner()
        r = designer.design(ss.dot_bracket(), ss.sequence())
        ss.replace_sequence(r[0].sequence)
        new_mg.replace_helix_sequence(ss)
        #new_mg.write_pdbs()

    def test_get_pose(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        new_mg = mg.copy()
        new_mg.replace_ideal_helices()
        ss = new_mg.designable_secondary_structure()

    def test_secondary_structure(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()
        for n in mt.tree.nodes:
            mg.add_motif(n.data)
        #mg.write_pdbs("org")
        mg.replace_ideal_helices()

        ss = mg.designable_secondary_structure()

        designer = sd.SequenceDesigner()
        r = designer.design(ss.dot_bracket(), ss.sequence())
        ss.replace_sequence(r[0].sequence)
        mg.replace_helix_sequence(ss)
        #mg.write_pdbs("new")
        #mg.merger.to_pdb("test.pdb")
        #mg.write_pdbs()

    def test_designable_secondary_structure(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        mg.replace_ideal_helices()
        #mg.write_pdbs()
        #for i, c in enumerate(mg.structure.chains()):
        #    c.to_pdb("c."+str(i)+".pdb")

        #mg.structure.to_pdb("test.pdb")
        mg.write_pdbs()
        ss = mg.designable_secondary_structure()

        designer = sd.SequenceDesigner()
        r = designer.design(ss.dot_bracket(), ss.sequence())
        ss.replace_sequence(r[0].sequence)
        mg.replace_helix_sequence(ss)
        #mg.write_pdbs("new")
        #mg.merger.to_pdb("test.pdb")

    def test_designable_secondary_structure_2(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        hairpin = rm.manager.mlibs['hairpin'].get_random()
        mg.add_motif(hairpin)

    def test_ss_real_case(self):
        rm.manager.add_motif("resources/motifs/tetraloop_receptor_min")
        mg = motif_graph.MotifGraph()
        mg.add_motif(rm.manager.get_motif(name="tetraloop_receptor_min",
                                          end_name="A228-A246"))
        mg.add_motif(rm.manager.get_motif(name="HELIX.IDEAL.20"), parent_end_name="A221-A252")
        mg.add_motif(rm.manager.get_motif(name="tetraloop_receptor_min",
                                          end_name="A221-A252"))
        mg.replace_ideal_helices()
        ss = mg.designable_secondary_structure()
        build.fill_basepairs_in_ss(ss)
        mg.replace_helix_sequence(ss)

    def test_to_tree(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)

        mt2 = motif_topology.graph_to_tree(mg)

    def test_from_ss(self):
        builder = build.BuildSecondaryStructure()
        ss_p = builder.build_helix(10)
        ssg  = secondary_structure_graph.graph_from_pose(ss_p)

    def test_topology_to_str(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mg = motif_graph.MotifGraph()

        for n in mt.tree.nodes:
            mg.add_motif(n.data)


        s = mg.topology_to_str_new()
        new_mg = motif_graph.MotifGraph(top_str=s)
        if len(mg) != len(new_mg):
            self.fail("did not get the correct number of nodes")

    def test_topology_to_str_2(self):
        #load motif graph with all atoms
        base_dir = settings.UNITTEST_PATH + "resources/motif_graph/"
        f = open(base_dir+"test.mg")
        l = f.readline()
        f.close()

        mg = motif_graph.MotifGraph(mg_str=l)

        #load motif graph from topology
        base_dir = settings.UNITTEST_PATH + "resources/motif_graph/"
        f = open(base_dir+"test.top")
        l = f.readline()
        f.close()

        mg2 = motif_graph.MotifGraph(top_str=l)

        atoms1 = mg.get_structure().structure.atoms()
        atoms2 = mg2.get_structure().structure.atoms()

        for i in range(len(atoms1)):
            if not numerical.are_atom_equal(atoms1[i], atoms2[i]):
                self.fail("atoms are not equal")

    def test_topology_to_str_3(self):
        #load motif graph with all atoms
        base_dir = settings.UNITTEST_PATH + "resources/motif_graph/"
        f = open(base_dir+"test_added.mg")
        l = f.readline()
        f.close()

        mg = motif_graph.MotifGraph(mg_str=l)

        #load motif graph from topology
        base_dir = settings.UNITTEST_PATH + "resources/motif_graph/"
        f = open(base_dir+"test_added.top")
        l = f.readline()
        f.close()

        mg2 = motif_graph.MotifGraph(top_str=l)

        atoms1 = mg.get_structure().structure.atoms()
        atoms2 = mg2.get_structure().structure.atoms()

        for i in range(len(atoms1)):
            if not numerical.are_atom_equal(atoms1[i], atoms2[i]):
                self.fail("atoms are not equal")

    def test_to_str(self):
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
            if not numerical.are_atom_equal(atoms1[i], atoms2[i]):
                self.fail("atoms are not equal")

    def test_to_str_2(self):
        builder = build.BuildMotifGraph()
        mg = builder.build(3)

        mg.get_node(0).index = 5
        mg.aligned[5] = 0
        mg.get_node(1).index = 7
        mg.aligned[7] = 1
        s = mg.to_str()
        mg_new = motif_graph.MotifGraph(mg_str=s)
        indices = []
        for n in mg.graph:
            indices.append(n.index)

        if " ".join(str(x) for x in indices) != "5 7 2":
            self.fail("to_str did not work")

    def test_replace_helix_sequence(self):
        builder = build.BuildMotifGraph()
        mg = builder.build(5)
        mg.replace_ideal_helices()

        dss = mg.designable_secondary_structure()
        secondary_structure_tools.fill_basepairs_in_ss(dss)
        mg.replace_helix_sequence(dss)

        mg_new = motif_graph.MotifGraph()
        for n in mg.graph:
            mg_new.add_motif(m_name=n.data.name, m_end_name=n.data.ends[0].name())

        rs1 = mg.get_structure()
        rs2 = mg_new.get_structure()

        if not is_equal.are_atoms_equal(rs1.structure.atoms(),
                                        rs2.structure.atoms()):
            self.fail("structures are not the same")

    def test_replace_helix_sequence_2(self):
        #load motif graph from topology
        base_dir = settings.UNITTEST_PATH + "resources/motif_graph/"
        f = open(base_dir+"mg_build.top")
        l = f.readline()
        f.close()

        mg = motif_graph.MotifGraph(top_str=l)
        mg.replace_ideal_helices()

        dss = mg.designable_secondary_structure()
        secondary_structure_tools.fill_basepairs_in_ss(dss)
        mg.replace_helix_sequence(dss)

        for n in graph.transverse_graph(mg.graph, 0):
            print n.index


def main():
    unittest.main()

if __name__ == '__main__':
    main()














