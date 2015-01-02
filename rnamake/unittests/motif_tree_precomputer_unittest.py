import unittest
import os
import rnamake.motif_tree_precomputer
import rnamake.resource_manager
import rnamake.motif_tree

class MotifTreePrecomputerUnittest(unittest.TestCase):

    def test_creation(self):
        mtp = rnamake.motif_tree_precomputer.MotifTreePrecomputer()

    def test_sqlite3_output(self):
        output = rnamake.motif_tree_precomputer.MotifTreePrecomputerSqlite3Output('test')
        os.remove("test_state.db")
        os.remove("test_data.db")

    def test_text_output(self):
        output = rnamake.motif_tree_precomputer.MotifTreePrecomputerTextOutput('test')
        os.remove("test.new.me")

    def test_precompute(self):
        rm = rnamake.resource_manager.ResourceManager()
        mtp = rnamake.motif_tree_precomputer.MotifTreePrecomputer()
        m = rm.get_motif("HELIX.IDEAL")
        mtp._precompute(m, 0, 0, 0, 0)
        try:
            os.remove("test.new.me")
        except:
            pass

    def test_precompute_motif(self):
        rm = rnamake.resource_manager.ResourceManager()
        m = rm.get_motif("TWOWAY.1GID.2")
        mtp = rnamake.motif_tree_precomputer.MotifTreePrecomputer()
        mtp.precompute_motif(m)


    def test_motif_orientation(self):
        mt = rnamake.motif_tree.MotifTree()
        mt.write_pdbs()
        rm = rnamake.resource_manager.ResourceManager()
        m = rm.get_motif("TWOWAY.1GID.2")
        hm = rm.get_motif("HELIX.IDEAL.9")
        motifs = []
        origin = mt.nodes[0].motif.ends[0].d()
        # print origin
        for end_index in (0, 1):
            for flip in (0, 1):
                for helix_end_index in (0, 1):
                    mt.add_motif(hm, end_index=helix_end_index)
                    mt.add_motif(m, end_index=end_index, end_flip=flip)
                    mt.add_motif(hm)
                    if len(mt.nodes) != 4:
                        #print helix_end_index, end_index, flip
                        mt.remove_node_level()
                        continue
                    motifs.append(mt.get_pose())
                    mt.remove_node_level()
        # for i, m in enumerate(motifs):
            # m.to_pdb("motif."+str(i)+".pdb")
           #  print m.ends[0].d(), m.ends[1].d()





def main():
    unittest.main()

if __name__ == '__main__':
    main()
