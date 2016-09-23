from rnamake.eternabot import sequence_designer
from rnamake import resource_manager as rm
from rnamake import motif_state_tree, monte_carlo
from rnamake import basepair

from collections import defaultdict
import random

class OptimizedSequence(object):
    def __init__(self):
        pass

class SequenceOptimizer(object):
    def __init__(self):
        pass

    def get_optimized_sequences(self, mg, node_i, node_j, end_i, end_j):
        dss = mg.designable_secondary_structure()
        org_dss = mg.designable_secondary_structure()
        node_name = mg.get_node(node_j).data.name



class SequenceOptimizer3D(object):
    def __init__(self):
        self.opt_nodes = []
        self.bps = {}

    class _DesignableBP(object):
        def __init__(self, bp):
            self.bp = bp
            self.static = 0
            self.state = ""
            self.last_state = ""
            self.end_0_node = -1
            self.end_1_node = -1

    def _update_designable_bp(self, d_bp, mst, new_bp_state):
        d_bp.state = new_bp_state

        if d_bp.end_0_node != -1:
            n = mst.get_node(d_bp.end_0_node)
            spl = n.data.cur_state.name.split("=")
            new_name = new_bp_state + "=" + spl[1]
            mst.replace_state(n.index, rm.manager.get_state(name=new_name))

        if d_bp.end_1_node != -1:
            n = mst.get_node(d_bp.end_1_node)
            spl = n.data.cur_state.name.split("=")
            new_name = spl[0] + "=" + new_bp_state
            mst.replace_state(n.index, rm.manager.get_state(name=new_name))


    def get_optimized_sequences(self, mt, target_bp, beads=None):
        mst = motif_state_tree.MotifStateTree(mt=mt)
        #mt.to_pdb("test.pdb")
        #mt.write_pdbs()

        for n in mt:
            if n.data.name != "HELIX.IDEAL":
                continue

            if n.parent is None:
                designable_bp_1 = self._DesignableBP(n.data.ends[0])
                designable_bp_1.end_0_node = n.index

                self.bps[n.data.ends[0]] = designable_bp_1

                designable_bp_2 = self._DesignableBP(n.data.ends[1])
                designable_bp_2.end_1_node = n.index

                self.bps[n.data.ends[1]] = designable_bp_2

            else:

                pei = n.parent_end_index()
                parent_end = n.parent.data.ends[pei]
                if parent_end in self.bps:
                    self.bps[parent_end].end_0_node = n.index
                else:
                    designable_bp_1 = self._DesignableBP(n.data.ends[0])
                    static_state = parent_end.res1.name + parent_end.res2.name
                    designable_bp_1.static = 1
                    designable_bp_1.state = static_state
                    designable_bp_1.end_0_node = n.index
                    self.bps[n.data.ends[0]] = designable_bp_1

                child = None
                if n.children[1] is not None:
                    if n.children[1].data.name != "HELIX.IDEAL":
                        child = n.children[1]


                if child:
                    designable_bp_2 = self._DesignableBP(n.data.ends[1])
                    child_end = child.data.ends[0]
                    static_state = child_end.res1.name + child_end.res2.name
                    designable_bp_2.static = 1
                    designable_bp_2.state = static_state
                    designable_bp_2.end_1_node = n.index
                    self.bps[n.data.ends[1]] = designable_bp_2
                else:
                    designable_bp = self._DesignableBP(n.data.ends[1])
                    designable_bp.end_1_node = n.index
                    self.bps[n.data.ends[1]] = designable_bp



        possible_bps = ["AU", "UA", "GC", "CG"]


        for bp in self.bps.values():
            if bp.static == 0:
                bp.state = random.choice(possible_bps)

        nodes_and_names = {}
        for bp in self.bps.values():

            if bp.end_0_node not in nodes_and_names:
                nodes_and_names[bp.end_0_node] = [bp.state, ""]
            else:
                nodes_and_names[bp.end_0_node][0] = bp.state

            if bp.end_1_node not in nodes_and_names:
                nodes_and_names[bp.end_1_node] = ["", bp.state]
            else:
                nodes_and_names[bp.end_1_node][1] = bp.state

        for index, names in nodes_and_names.iteritems():
            if index == -1:
                continue
            name = "=".join(names)
            mst.replace_state(index, rm.manager.get_state(name=name))

        designable_bps = []
        for bp in self.bps.values():
            if bp.static == 0:
                designable_bps.append(bp)

        target_state = target_bp.state()
        target_state_flip = target_bp.state()
        target_state_flip.flip()

        last_score = target_state.diff(mst.last_node().data.cur_state.end_states[1])

        mc = monte_carlo.MonteCarlo()
        mc.temperature = 0.5
        best = 1000
        best_states = []
        for i in range(500):
            d_bp = random.choice(designable_bps)

            d_bp.last_state = d_bp.state
            new_bp_state = random.choice(possible_bps)

            self._update_designable_bp(d_bp, mst, new_bp_state)

            new_score = target_state.diff(mst.last_node().data.cur_state.end_states[1])


            if mc.accept(last_score, new_score):
                last_score = new_score
            else:
                self._update_designable_bp(d_bp, mst, d_bp.last_state)


            if best > new_score:
                best = new_score
                best_states = []
                for n in mst:
                    best_states.append([n.data.cur_state.name,
                                        n.data.cur_state.end_names[0]])



        return best, best_states
        #mst.to_motif_tree().to_pdb("opt.pdb", renumber=1, close_chain=1)
        #print mst.to_motif_tree().secondary_structure().sequence()


