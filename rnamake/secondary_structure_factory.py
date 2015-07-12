import ss_tree
import secondary_structure

class StructureSecondaryFactory(object):
    def __init__(self):
        pass

    def _get_basepairs(self, sstree, ss):
        basepairs, ends = [], []
        for n in sstree:
            if n.data.type == ss_tree.SS_Type.SS_BP:
                res1 = n.data.ss_chains[0].residues[0]
                res2 = n.data.ss_chains[1].residues[0]

                bp = secondary_structure.Basepair(res1, res2)
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
        return basepairs, ends

    def _get_elements(self, sstree, ss):
        elements = { 'ALL' : [] }

        for n in sstree:
            if n.parent_index() == -1 or n.parent.data.type == ss_tree.SS_Type.SS_SEQ_BREAK:
                continue

            if n.data.type == ss_tree.SS_Type.SS_BP and \
               n.parent.data.type == ss_tree.SS_Type.SS_BP:
                nodes = [n,n.parent]
                bp_nodes = [n,n.parent]
                type_name = "BP_STEP"

            elif n.data.type != ss_tree.SS_Type.SS_SEQ_BREAK and \
                n.parent.data.type == ss_tree.SS_Type.SS_BP:
                nodes = [n.parent, n]
                bp_nodes = [n.parent]
                type_name = n.data.what()[3:]
                for c in n.children:
                    if c is None:
                        continue
                    if c.data.type == ss_tree.SS_Type.SS_BP or \
                       c.data.type == ss_tree.SS_Type.SS_PSEUDO_BP:
                        nodes.append(c)
                        bp_nodes.append(c)
                    else:
                        print n.data.what(), c.data.what()
                        raise ValueError("unexpected connectivity in ss_tree")
            else:
                continue


            chains = self._get_chains(ss, nodes)
            ends = []
            for bp_n in bp_nodes:
                bp = ss.get_bp(bp_n.data.ss_chains[0].residues[0],
                               bp_n.data.ss_chains[1].residues[0])
                ends.append(bp)
            e = secondary_structure.SecondaryStructureMotif(type_name, ends, chains)
            bps = []
            res = {r : 1 for r in e.residues() }
            for bp in ss.basepairs:
                if bp.res1 in res and bp.res2 in res:
                    bps.append(bp)

            e.basepairs = bps
            if type_name not in elements:
                elements[type_name] = []

            elements[type_name].append(e)
            elements['ALL'].append(e)
        return elements

    def _get_chains(self, ss, nodes):
        res = []
        for n in nodes:
            for c in n.data.ss_chains:
                res.extend(c.residues)
        res.sort(key=lambda x: x.num)
        chains = []
        c_res = []
        last = -1
        for r in res:
            is_chain_start = 0
            for c in ss.chains:
                if c.first() == r:
                    is_chain_start = 1
                    break

            if last == -1:
                pass
            elif last+1 != r.num or is_chain_start:
                chains.append(secondary_structure.Chain(c_res))
                c_res = []
            c_res.append(r)
            last = r.num
        if len(c_res) > 0:
            chains.append(secondary_structure.Chain(c_res))
        return chains

    def get_structure(self, sequence=None, dot_bracket=None, base_ss=None):
        if   sequence is not None and dot_bracket is not None:
            sstree = ss_tree.SS_Tree(sequence, dot_bracket)
        elif base_ss is not None:
            sstree = ss_tree.SS_Tree(ss=base_ss)
        else:
            raise ValueError("supply sequence and dot_bracket strings or a" + \
                             " SecondaryStructure object")
        ss = sstree.ss

        basepairs, ends = self._get_basepairs(sstree, ss)
        ss.ends = ends
        ss.basepairs = basepairs
        ss.elements = self._get_elements(sstree, ss)

        return ss

factory = StructureSecondaryFactory()

def ss_id_to_ss_tree(ss_id):
    ss = ""
    seq = ""
    spl = ss_id.split("_")

    for i in range(0, len(spl)-1, 2):
        seq += spl[i]
        for e in spl[i+1]:
            if   e == "L":
                ss += "("
            elif e == "R":
                ss += ")"
            elif e == "U":
                ss += "."
            else:
                raise ValueError("unexpected symbol in ss_id")

        if i != len(spl)-2:
            seq += "+"
            ss += "+"

    return ss_tree.SS_Tree(ss, seq)

