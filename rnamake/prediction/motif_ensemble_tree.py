import rnamake.motif_tree_state as motif_tree_state
import rnamake.motif_tree_state_tree as motif_tree_state_tree
import rnamake.motif_type as motif_type
import rnamake.util as util
import rnamake.resource_manager as resource_manager
import motif_ensemble
import random
import copy

class MotifEnsembleTree(object):
    def __init__(self, mt=None, ensemble=None):
        if ensemble is None:
            head = self._default_head()
        else:
            head = MotifEnsembleTreeNode(ensemble, None, 0)
        self.nodes = [ head ]
        self.last_node = head

        if mt is not None:
            self._setup_from_mt(mt)

    def _setup_from_mt(self, mt):
        for i, n in enumerate(mt.nodes):
            if i == 0:
                continue

            parent_index = n.parent().index
            parent_end_index = n.parent_end_index()
            end_index = n.motif_end_index(n.parent())

            if n.motif.name[3] == "=":
                me = motif_ensemble.MotifEnsemble(n.motif.name)
            else:
                mts_name = n.motif.name + "-" + str(end_index)
                mts = resource_manager.manager.get_state(mts_name)
                me = mts_to_me(mts)

            node = self.add_ensemble(me, self.nodes[parent_index],
                                     parent_end_index)

            if node is None:
                raise ValueError("could not translate motiftree to motifensembletree")


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
            mtst = motif_tree_state_tree.MotifTreeStateTree()
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


class MTSTtoMETConverter(object):
    def __init__(self):
        pass

    def _get_start_chain(self, n):
        residues = n.motif.residues()
        closest_to_5prime = 10000
        start_chain = None
        res = None

        for r in residues:
            r_new = self.p.get_residue(uuid=r.uuid)
            if r_new is None:
                continue

            for i, c in enumerate(self.p.chains()):
                if r_new in c.residues:
                    index = c.residues.index(r_new)
                    if index < closest_to_5prime:
                        closest_to_5prime = index
                        res = r_new
                        start_chain = c

        return start_chain

    def _get_chain_pos(self, r):
        for i, c in enumerate(self.p.chains()):
            if r in c.residues:
                return i

    def _add_helix(self, n, chain, motif_bp, start_index, flip):
        bps_str = []
        if motif_bp is not None:
            bps_str.append(motif_bp)
        chain_pos = self._get_chain_pos(chain.first())

        start_pos = 10000
        for r in n.motif.residues():
            r_new = self.p.get_residue(uuid=r.uuid)
            if r_new is None:
                continue
            if r_new in chain.residues:
                index = chain.residues.index(r_new)
                if index < start_pos:
                    start_pos = index

        for i, r in enumerate(chain.residues):
            if i < start_pos:
                continue
            r_new = n.motif.get_residue(uuid=r.uuid)
            if r_new is None:
                break
            res1 = self.dseq[r.num-1 + chain_pos]
            bps = self.p.get_basepair(uuid1=r.uuid)
            for bp in bps:
                if bp.res1 == r:
                    res2 =  self.dseq[bp.res2.num-1 + self._get_chain_pos(bp.res2)]
                    bps_str.append(res1+res2)
                else:
                    res2 =  self.dseq[bp.res1.num-1 + self._get_chain_pos(bp.res1)]
                    bps_str.append(res2+res1)
        steps = []
        for i in range(1, len(bps_str)):
            steps.append(bps_str[i-1]+"="+bps_str[i])

        index = 0
        for c in n.connections:
            partner = c.partner(n)
            if partner.index < n.index:
                index = partner.motif.ends.index(c.motif_end(partner))
                break
        me = motif_ensemble.MotifEnsemble(steps[0], start_index, flip)
        self.met.add_ensemble(me, parent_end_index=index)

        for s in steps[1:]:
            me = motif_ensemble.MotifEnsemble(s, start_index, flip)
            self.met.add_ensemble(me)

        return bps_str[-1]

    def _get_next_bp(self, n, next_node, chain, flip=0):
        if not flip:
            bp = n.connection(next_node).motif_end(next_node)
        else:
            bp = n.connection(next_node).motif_end(n)
        first_m_bp = None
        res1 = None
        best_index = 10000
        for r in bp.residues():
            r_new = self.p.get_residue(uuid=r.uuid)
            if r_new in chain.residues:
                index = chain.residues.index(r_new)
                if index < best_index:
                    res1 = r
                    best_index = index

        if res1 == bp.res1:
            first_m_bp = bp.res1.rtype.name[0] +  bp.res2.rtype.name[0]
        else:
            first_m_bp = bp.res2.rtype.name[0] +  bp.res1.rtype.name[0]
        return first_m_bp

    def convert(self, mtst, start_pos=1, pose=None, seq=None, ss=None, debug=0):
        if not pose:
            self.p = mtst.to_pose()
        else:
            self.p = pose
        mts_lib = motif_tree_state.MotifTreeStateLibrary(motif_type.TWOWAY)
        if debug:
            self.dseq = self.p.sequence()
            score = 0
        else:
            self.dseq, score = self.p.optimized_sequence(ss)
        print self.dseq

        start_chain = self._get_start_chain(self.p.nodes[start_pos])
        self.met = MotifEnsembleTree()

        for i, n in enumerate(self.p.nodes):
            if i < start_pos:
                continue
            #ne = motif_tree_state.parse_db_name(mtst.nodes[i].mts.name)

            motif_bp = None
            if n.motif.mtype == motif_type.HELIX and "TWOWAY" not in n.motif.name:
                if i > start_pos:
                    motif_bp = self._get_next_bp(self.p.nodes[i-1], n, start_chain, 1)

                flip = mtst.nodes[i].mts.flip
                start_index = mtst.nodes[i].mts.start_index
                last_bp = self._add_helix(n, start_chain, motif_bp,
                                          start_index,
                                          flip)

                if i < len(self.p.nodes)-1:
                    next_bp = self._get_next_bp(n, self.p.nodes[i+1], start_chain)
                    last_step = last_bp + "=" + next_bp
                    me = motif_ensemble.MotifEnsemble(last_step,
                                                      start_index,
                                                      flip)
                    self.met.add_ensemble(me)

                elif len(n.connections) > 1:
                    other_node = None
                    for c in n.connections:
                        partner = c.partner(n)
                        if partner == self.p.nodes[i-1]:
                            continue
                        other_node = partner
                    next_bp = self._get_next_bp(n, other_node, start_chain)
                    last_step = last_bp + "=" + next_bp
                    me = motif_ensemble.MotifEnsemble(last_step,
                                                      start_index,
                                                      flip)
                    self.met.add_ensemble(me)



            else:
                me =  mts_to_me(mtst.nodes[i].mts)
                self.met.add_ensemble(me)
                #flipped = self._add_mts(mtst, i, mts_lib)

        #mtst.to_pdb('test.pdb')
        mtst2 = self.met.get_mtst()
        #mtst2.to_pdb("test2.pdb")
        return mtst2, score

def mts_to_me(mts):
    me = motif_ensemble.MotifEnsemble()
    ms = motif_ensemble.MotifState(mts, 1.0)
    me.motif_states.append(ms)
    return me

