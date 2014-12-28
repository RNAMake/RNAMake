import unittest
import rnamake.motif_tree_merger
import rnamake.motif_tree
import rnamake.resource_manager
import rnamake.settings
import rnamake.motif_tree_state
import instance

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
        # print len(mt.nodes)
        mt.write_pdbs()
        mt._find_other_connections_to_head()
        merger = rnamake.motif_tree_merger.MotifTreeMerger()
        pose = merger.merge(mt,include_head=1)
        pose.to_pdb()




def main():
    unittest.main()

if __name__ == '__main__':
    main()
