import unittest
import rnamake.motif_tree as motif_tree
import rnamake.motif_factory as motif_factory
import rnamake.motif_state_tree as motif_state_tree
import rnamake.resource_manager as rm
import build

class MotifStateTreeUnittest(unittest.TestCase):

    def _test_creation(self):
        builder = build.BuildMotifTree()
        mt = builder.build(10)

        mst = motif_state_tree.MotifStateTree()

        for i, n in enumerate(mt):
            if i == 0:
                mst.add_state(rm.manager.get_state(n.data.name))
            else:
                parent_index = 1000
                parent_end_index = -1
                for c in n.connections:
                    if c is None:
                        continue
                    if c.partner(n.index).index < parent_index:
                        parent_index =  c.partner(n.index).index
                        parent_end_index = c.end_index(c.partner(n.index).index)

                if parent_index == 1000:
                    raise ValueError("did not convert motif tree to motif state tree properly")
                mst.add_state(rm.manager.get_state(n.data.name), parent_index, parent_end_index)

    def test_creation_from_mt(self):
        builder = build.BuildMotifTree()
        mt = builder.build(10)
        mst = motif_state_tree.MotifStateTree(mt)

    def test_to_mt(self):
        builder = build.BuildMotifTree()
        mt = builder.build(10)
        mt.to_pdb("test1.pdb")
        mst = motif_state_tree.MotifStateTree(mt)

        mt2 = mst.to_motif_tree()
        mt2.to_pdb("test2.pdb")



def main():
    unittest.main()

if __name__ == '__main__':
    main()