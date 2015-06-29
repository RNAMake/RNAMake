import unittest
import build
import rnamake.sqlite_library as sqlite_library
import rnamake.secondary_structure as secondary_structure
import rnamake.resource_manager as rm
import rnamake.motif_tree as motif_tree
import rnamake.motif_factory as motif_factory
import rnamake.ss_tree as ss_tree

class SecondaryStructureUnittest(unittest.TestCase):

    def test_assign_secondary_structure(self):
        builder = build.BuildMotifTree()
        mt = builder.build()
        #for n in mt:
        #    print n.data.name

        p = mt.to_pose()
        #print p.secondary_structure()
        #print p.sequence()
        #mt.write_pdbs()

    def test_parse(self):
        m = rm.manager.get_motif("NWAY.1S72.18-01551-01634")
        ss = m.secondary_structure()
        seq = m.sequence()

        sstree = ss_tree.SS_Tree(ss, seq)

        chains = []
        basepairs = []
        ends = []
        all_res = []
        residues = []

        for i in range(len(seq)):
            if seq[i] != "&":
                r = secondary_structure.Residue(seq[i], i)
                residues.append(r)
                all_res.append(r)
            else:
                chains.append(secondary_structure.Chain(residues))
                residues = []

        if len(residues) > 0:
            chains.append(secondary_structure.Chain(residues))


        for n in sstree:
            if n.data.type == ss_tree.SS_Type.SS_BP:
                res1_i = n.data.ss_data[0].bounds[0]
                res2_i = n.data.ss_data[1].bounds[1]
                bp_res = []
                for r in all_res:
                    if r.num == res1_i or r.num == res2_i:
                        bp_res.append(r)
                if len(bp_res) != 2:
                    raise ValueError("incorrect number of resiudes in ss_bp")

                bp = secondary_structure.Basepair(bp_res[0], bp_res[1])
                basepairs.append(bp)
                if n.parent is None:
                    ends.append(bp)
                    continue
                children = []
                for c in n.children:
                    if c is None:
                        continue
                    children.append(c)
                if len(children) == 1:
                    if children[0].data.type == ss_tree.SS_Type.SS_SEQ_BREAK:
                        ends.append(bp)


        struct = secondary_structure.Structure(chains, basepairs)
        struct.ends = ends

        print ends[0].res1.num, ends[0].res2.num


        print seq
        print ss

        ss_chains = struct.reorient_ss_and_seq(ends[2].res1.num, ends[2].res2.num)

        print struct.sequence()
        print struct.secondary_structure()



def main():
    unittest.main()

if __name__ == '__main__':
    main()