import unittest
import rnamake.motif_tree_merger
import rnamake.motif_tree
import rnamake.resource_manager
import rnamake.settings
import rnamake.motif_tree_state
import instance
import numerical
import redesign.chain_closure

class MotifTreeMergerUnittest(unittest.TestCase):

    def test_creation(self):
        mt = instance.simple_mt()
        merger = rnamake.motif_tree_merger.MotifTreeMerger()
        pose = merger.merge(mt)
        s1 = "C1 C2 C3 "
        s2 = "G4 G5 G6 "
        if pose.chains()[0].list_res() != s1:
            self.fail()
        if pose.chains()[1].list_res() != s2:
            self.fail()

    def test_creation_2(self):
        mt = instance.simple_mt_with_head()
        merger = rnamake.motif_tree_merger.MotifTreeMerger()
        pose = merger.merge(mt,include_head=1)

    def test_merge(self):
        mt = instance.simple_mt_helix(10)
        merger = rnamake.motif_tree_merger.MotifTreeMerger()
        pose = merger.merge(mt)

    def test_merge_2(self):
        mt =  instance.simple_mt_with_head()
        rm = rnamake.resource_manager.ResourceManager()
        m = rm.get_motif("HELIX.IDEAL")
        mt.add_motif(m, parent=mt.nodes[0])
        merger = rnamake.motif_tree_merger.MotifTreeMerger()
        pose = merger.merge(mt,include_head=0)

    def test_connection(self):
        #TODO write function to see everything on a chain is in order
        rm = rnamake.resource_manager.ResourceManager()
        path = rnamake.settings.UNITTEST_PATH + "/resources/motifs"
        rm.add_lib_path(path)
        s = "HELIX.LE.16-0-0-0-0-1-1,TWOWAY.1GID.2-0-0-0-0-1-0,"+\
        "HELIX.LE.11-1-0-1-0-0-1,TWOWAY.2HW8.0-1-0-1-0-0-0,HELIX.LE.7-1-0-1-0-0-0,"+\
        "TWOWAY.2HW8.0-1-0-1-0-0-0"
        mt = rnamake.motif_tree.MotifTree(rm.get_motif("tetraloop_receptor_min"))
        parent_index=0
        for i, db_name in enumerate(s.split(",")):
            ne = rnamake.motif_tree_state.parse_db_name(db_name)
            m = rm.get_motif(ne.motif_name)
            mt.add_motif(m, end_index=ne.start_index, end_flip=ne.flip_direction,
                             parent_index=parent_index)
            parent_index = ne.end_index
        # mt.write_pdbs()
        mt._find_other_connections_to_head()
        merger = rnamake.motif_tree_merger.MotifTreeMerger()
        pose = merger.merge(mt,include_head=1)
        # pose.to_pdb()

    def test_create_coord_system(self):
        rm = rnamake.resource_manager.ResourceManager()
        m = rm.get_motif("HELIX.IDEAL")
        res1 = m.residues()[0]
        atoms = [res1.get_atom("O3'"),res1.get_atom("C3'"),res1.get_atom("C4'")]
        matrix = rnamake.motif_tree_merger.create_coord_system(atoms)
        m_old = redesign.chain_closure.create_coord_system(atoms)
        if not numerical.are_matrices_equal(matrix, m_old):
            self.fail()

    def test_get_virtal_atom(self):
        rm = rnamake.resource_manager.ResourceManager()
        m = rm.get_motif("HELIX.IDEAL")
        res1 = m.residues()[0]
        atoms = [res1.get_atom("O3'"),res1.get_atom("C3'"),res1.get_atom("C4'")]
        ovl1 =  rnamake.motif_tree_merger.virtual_atom("OVL1", 1.606497,
                                                        60.314519, 0.0,
                                                        atoms)
        old_ovl1 = redesign.chain_closure.make_virtual_atom("OVL1", 1.606497,
                                                           60.314519, 0.0,
                                                           atoms)

        if not numerical.are_atom_equal(ovl1, old_ovl1):
            self.fail()

    def test_chain_closure(self):
        rm = rnamake.resource_manager.ResourceManager()
        path = rnamake.settings.UNITTEST_PATH + "/resources/motifs"
        rm.add_lib_path(path)
        s = "HELIX.LE.16-0-0-0-0-1-1,TWOWAY.1GID.2-0-0-0-0-1-0,"+\
        "HELIX.LE.11-1-0-1-0-0-1,TWOWAY.2HW8.0-1-0-1-0-0-0,HELIX.LE.7-1-0-1-0-0-0,"+\
        "TWOWAY.2HW8.0-1-0-1-0-0-0"
        mt = rnamake.motif_tree.MotifTree(rm.get_motif("tetraloop_receptor_min"))
        parent_index=0
        for i, db_name in enumerate(s.split(",")):
            ne = rnamake.motif_tree_state.parse_db_name(db_name)
            m = rm.get_motif(ne.motif_name)
            mt.add_motif(m, end_index=ne.start_index, end_flip=ne.flip_direction,
                             parent_index=parent_index)
            parent_index = ne.end_index
        # mt.write_pdbs()
        mt._find_other_connections_to_head()
        merger = rnamake.motif_tree_merger.MotifTreeMerger()
        pose = merger.merge(mt,include_head=1)

        for c in pose.chains():
            rnamake.motif_tree_merger.close_chain(c)
        pose.to_pdb()







def main():
    unittest.main()

if __name__ == '__main__':
    main()
