import rnamake.motif_tree_state as motif_tree_state
import rnamake.motif_type as motif_type
import motif_ensemble
import random

class MotifEnsembleTree(object):
    def __init__(self, ensemble=None):
        if ensemble is None:
            head = self._default_head()
        else:
            head = MotifEnsembleTreeNode(ensemble, None, 0)
        self.nodes = [ head ]
        self.last_node = head

    def add_ensemble(self, ensemble, parent=None, parent_end_index=None):
        if parent is None:
            parent = self.last_node
        if parent_end_index is None:
            parent_end_index = parent.available_ends()[0]
        node = MotifEnsembleTreeNode(ensemble, parent, len(self.nodes))
        parent.add_child(node, parent_end_index)
        self.nodes.append(node)
        self.last_node = node
        return node

    def _default_head(self):
        me = motif_ensemble.MotifEnsemble()
        mts = motif_tree_state.ref_mts()
        ms = motif_ensemble.MotifState(mts, 1.0)
        me.motif_states.append(ms)
        return MotifEnsembleTreeNode(me, None, 0)

    def get_mtst(self):
        mtst = None
        mts_name = self.nodes[0].motif_ensemble.motif_states[0].mts.name
        if mts_name == "start":
            mtst = motif_tree_state.MotifTreeStateTree(sterics=0)
        else:
            mts =  self.nodes[0].motif_ensemble.motif_states[0].mts
            mtst = motif_tree_state.MotifTreeStateTree(mts)
        open_nodes = [ self.nodes[0] ]
        while len(open_nodes) != 0:
            current = open_nodes.pop(0)
            for i, n in enumerate(current.children):
                if n is None:
                    continue
                mts = n.motif_ensemble.motif_states[0].mts
                parent = mtst.nodes [ current.index ]
                parent_end = parent.states[ i ]
                mstn = mtst.add_state(mts, parent, parent_end)
                if mstn is None:
                    raise ValueError("cannot build mstn based on motif ensemble \
                                     tree topology")
                open_nodes.append(n)
        return mtst


class MotifEnsembleTreeNode(object):
    def __init__(self, motif_ensemble, parent, index):
        self.motif_ensemble, self.parent = motif_ensemble, parent
        self.index = index
        self.children = [None for e in self.motif_ensemble.motif_states[0].mts.end_states]

    def available_ends(self):
        ends = []
        for i, c in enumerate(self.children):
            if c is None:
                ends.append(i)
        return ends

    def add_child(self, node, end_direction):
        if self.children[end_direction] is not None:
            raise ValueError("cannot add child already is using this end")
        self.children[end_direction] = node

    def get_random_state(self):
        return random.choice(self.motif_ensemble.motif_states)



def _get_start_res(n, p):
    residues = n.motif.residues()
    closest_to_5prime = 10000
    start_chain = None
    chain_pose = 0
    res = None

    for r in residues:
        r_new = p.get_residue(uuid=r.uuid)
        if r_new is None:
            continue

        for i, c in enumerate(p.chains()):
            if r_new in c.residues:
                index = c.residues.index(r_new)
                if index < closest_to_5prime:
                    closest_to_5prime = index
                    res = r_new
                    start_chain = c
                    chain_pos = i

    return res, start_chain


def get_chain_pos(p, r):
    for i, c in enumerate(p.chains()):
        if r in c.residues:
            return i


def _add_helix(n, p, seq, chain, met, end, flip):
    bps_str = []
    chain_pos = 0
    for i, c in enumerate(p.chains()):
        if c == chain:
            chain_pos = i
            break

    for i, r in enumerate(chain.residues):
        r_new = n.motif.get_residue(uuid=r.uuid)
        if r_new is None:
            break
        res1 = seq[r.num-1 + chain_pos]
        bps = p.get_basepair(uuid1=r.uuid)
        for bp in bps:
            if bp.res1 == r:
                res2 =  seq[bp.res2.num-1 + get_chain_pos(p, bp.res2)]
            else:
                res2 =  seq[bp.res1.num-1 + get_chain_pos(p, bp.res1)]
        bps_str.append(res1+res2)

    steps = []
    for i in range(1, len(bps_str)):
        steps.append(bps_str[i-1]+"="+bps_str[i])


    print steps[0]
    me = motif_ensemble.MotifEnsemble(steps[0], end, flip)
    met.add_ensemble(me)
    for s in steps[1:]:
        break
        me = motif_ensemble.MotifEnsemble(s, end, 0)
        met.add_ensemble(me)

    return bps_str[-1]


def mts_to_me(mts):
    me = motif_ensemble.MotifEnsemble()
    ms = motif_ensemble.MotifState(mts, 1.0)
    me.motif_states.append(ms)
    return me


def mtst_to_met(mtst, start_pos=1):
    mtst.nodes_to_pdbs()
    print len(mtst.nodes)
    exit()
    p = mtst.to_pose()
    name_elements = motif_tree_state.parse_db_name(mtst.nodes[start_pos].mts.name)
    flip = 0
    start_index = 0
    print name_elements.start_index, name_elements.flip_direction
    if name_elements.start_index != 0 or name_elements.flip_direction != 0:
        flip = 1
    if  name_elements.start_index != 0 and name_elements.flip_direction != 0:
        flip = 0
        start_index = 1
    p.to_pdb('test.pdb')
    dseq = p.sequence()
    #dseq = p.optimized_sequence()
    print dseq

    res, start_chain = _get_start_res(p.nodes[start_pos], p)
    met =  MotifEnsembleTree()

    for i, n in enumerate(p.nodes):
        if i < start_pos:
            continue
        if n.motif.mtype == motif_type.HELIX:
            last_bp = _add_helix(n, p, dseq, start_chain, met,
                                 name_elements.start_index, name_elements.flip_direction)
            bp = n.connection(p.nodes[2]).motif_end(p.nodes[2])
            first_m_bp = None
            res1 = None
            for r in bp.residues():
                r_new = p.get_residue(uuid=r.uuid)
                if r_new in start_chain.residues:
                    res1 = r
            if res1 == bp.res1:
                first_m_bp = bp.res1.rtype.name[0] +  bp.res2.rtype.name[0]
            else:
                first_m_bp = bp.res2.rtype.name[0] +  bp.res1.rtype.name[0]
            last_step = last_bp + "=" + first_m_bp
            print last_step
            #me = motif_ensemble.MotifEnsemble(last_step,  name_elements.start_index, 0)
            #me = motif_ensemble.MotifEnsemble(last_step, 0, 0)
            #met.add_ensemble(me)
            break

        else:
            me = mts_to_me(mtst.nodes[i].mts)
            met.add_ensemble(me)

    for s in met.nodes[-2].motif_ensemble.motif_states[0].mts.end_states:
        if s is not None:
            print s.r
    for s in met.nodes[-1].motif_ensemble.motif_states[0].mts.end_states:
        if s is not None:
            print s.r


    mtst2 = met.get_mtst()
    mtst2.to_pdb("test2.pdb")




