import unittest
import random
import rnamake.motif_tree_state_path_optimizer
import rnamake.motif_tree_state as motif_tree_state
import rnamake.motif_type as motif_type
import rnamake.settings as settings
import rnamake.build_task_astar as build_task_astar
import rnamake.secondary_structure_tree as secondary_structure_tree
import rnamake.motif as motif
import rnamake.prediction.motif_ensemble_tree as motif_ensemble_tree
import rnamake.prediction.motif_ensemble as motif_ensemble
import rnamake.eternabot.sequence_designer as sequence_designer
import rnamake.util as util



def get_twoway_helix_mts_tree(size=4):
    twoways = motif_tree_state.MotifTreeStateLibrary(motif_type.TWOWAY)
    helixs = motif_tree_state.MotifTreeStateLibrary(motif_type.HELIX)
    me_libs = [helixs, twoways]
    mtst = motif_tree_state.MotifTreeStateTree()
    pos = 0
    i = 0
    count = 0
    while i < size:
        if i % 2 == 0:
            pos = 0
        else:
            pos = 1
        mts = random.choice(me_libs[pos].motif_tree_states)
        node = mtst.add_state(mts)
        if node is not None:
            i += 1
        count += 1
        if count > 1000:
            break
    return mtst

def mts_to_me(mts):
    me = motif_ensemble.MotifEnsemble()
    ms = motif_ensemble.MotifState(mts, 1.0)
    me.motif_states.append(ms)
    return me


def get_wc_flow():
    seq = 'CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG'
    ss  = '((((((....((((((((((((....))))))))))))....))))))'

    return seq, ss


def extract_steps_from_ss_tree(ss_tree):
    node = None
    for n in ss_tree.nodes:
        if n.ss_type == "Bulge":
            node = n
            break

    steps = []
    while node != None:
        if len(node.children) == 0:
            break
        c = node.children[0]
        if c.ss_type != "Basepair":
            break
        step = list(c.bp_type)
        if step[0] == "T":
            step[0] = "U"
        if step[1] == "T":
            step[1] = "U"
        steps.append("".join(step))
        node = c
    return steps


def get_step_motifs(steps):
    motif_names = []
    for i in range(1,len(steps)):
        full_step = steps[i-1] + "=" + steps[i]
        motif_names.append(full_step)
    return motif_names


def get_met(flow_ss_tree, chip_ss_tree):
    f = open('/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/simulate_tectos/tetraloop.str')
    lines = f.readlines()
    f.close()

    ggaa_motif = motif.str_to_motif(lines[0])
    gaaa_motif = motif.str_to_motif(lines[1])

    ggaa_state = motif_tree_state.motif_to_state(ggaa_motif, end_index=1,
                                                end_flip=1)
    gaaa_state = motif_tree_state.motif_to_state(gaaa_motif, end_index=0,
                                                end_flip=1)

    flow_steps = extract_steps_from_ss_tree(flow_ss_tree)

    motif_names = get_step_motifs(flow_steps[1:])

    met = motif_ensemble_tree.MotifEnsembleTree()
    #flow
    met.add_ensemble(motif_ensemble.MotifEnsemble("AU=GC", 0, 0))
    met.add_ensemble(mts_to_me(ggaa_state))
    met.add_ensemble(motif_ensemble.MotifEnsemble(motif_names[0], 0, 1), parent_end_index=2)
    for i in range(1, len(motif_names)):
        met.add_ensemble(motif_ensemble.MotifEnsemble(motif_names[i], 0, 1))
    #chip
    met.add_ensemble(mts_to_me(gaaa_state))

    return met


class MotifTreeStatePathOptimizerUnittest(unittest.TestCase):

    def get_test_case(self):
        problem = get_twoway_helix_mts_tree(5)
        path = settings.UNITTEST_PATH + "/resources/path_optimizer_prob.dat"
        f = open(path, 'w')
        for n in problem.nodes[1:]:
            f.write(n.mts.name + " ")
        f.write("\n")

        problem.to_pdb('prob.pdb')

        mtss = build_task_astar.MotifTreeStateSearch(max_solutions=1)
        solutions = mtss.search(problem.nodes[0].active_states()[0],
                    problem.last_node.active_states()[0])

        for i, s in enumerate(solutions):
            for n in s.path[1:]:
                 f.write(n.mts.name + " ")
            f.write("\n")
            s.to_mtst().to_pdb('solution.'+str(i)+'.pdb')
        f.close()


    def test(self):
        return
        #self.get_test_case()
        path = settings.UNITTEST_PATH + "/resources/path_optimizer_prob.dat"
        f = open(path)
        lines = f.readlines()
        f.close()
        twoways = motif_tree_state.MotifTreeStateLibrary(motif_type.TWOWAY)
        helixs = motif_tree_state.MotifTreeStateLibrary(motif_type.HELIX)
        mts_libs = [helixs, twoways]

        rm = resource_manager.ResourceManager()
        mtst = motif_tree_state.MotifTreeStateTree()
        all_mts = lines[0].split()
        pos = 0
        for i, mts in enumerate(all_mts):
            if i % 2 == 0:
                pos = 0
            else:
                pos = 1

            mtst.add_state(mts_libs[pos].get_state(mts))

        path = []
        all_mts = lines[1].split()
        mtst2 = motif_tree_state.MotifTreeStateTree()
        for i, mts in enumerate(all_mts):
            if i % 2 == 0:
                pos = 0
            else:
                pos = 1
            mtst2.add_state(mts_libs[pos].get_state(mts))

        mtspo = rnamake.motif_tree_state_path_optimizer.MotifTreeStatePathOptimizer()
        mtspo.optimize(mtst2)


    def test_tecto(self):
        return
        seq, ss = get_wc_flow()
        flow_ss_tree = secondary_structure_tree.SecondaryStructureTree(ss, seq)
        chip_ss_tree = secondary_structure_tree.SecondaryStructureTree(ss, seq)

        met = get_met(flow_ss_tree, chip_ss_tree)
        mtst = met.get_mtst()
        p = mtst.to_pose()
        chains = p.chains()
        chains.pop()
        basepairs = []
        res = []
        for c in chains:
            res.extend(c.residues)
        for bp in p.basepairs:
            if bp.res1 in res and bp.res2 in res:
                basepairs.append(bp)

        m = motif.Motif()
        m.structure._build_chains(res)
        m.structure._cache_coords()
        m.basepairs = basepairs
        m._cache_basepair_frames()
        m.setup_basepair_ends()

        mts = motif_tree_state.motif_to_state(m)
        mtst2 = motif_tree_state.MotifTreeStateTree(mts)

        twoways = motif_tree_state.MotifTreeStateLibrary(motif_type.TWOWAY)
        helixs = motif_tree_state.MotifTreeStateLibrary(motif_type.HELIX)
        mts_libs = [helixs, twoways]

        f = open('/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/motif_tree_path_refiner/out')
        lines = f.readlines()
        f.close()
        lines.pop(0)

        prime5_seq = 'CTAGGATATGG'
        prime3_seq = 'CCTAAGTCCTAG'
        prime5_ss  = '(((((((..(('
        prime3_ss  = '))...)))))))'

        #mtspo = rnamake.motif_tree_state_path_optimizer.MotifTreeStatePathOptimizer()
        designer = sequence_designer.SequenceDesigner()

        f = open('sequences.txt', 'w')

        pos = 0
        count = -1
        for l in lines:
            count += 1
            spl = l.split('|')
            mts_spl = spl[1].split()
            mts_spl.pop(0)
            for i, mts_name in enumerate(mts_spl):
                if i % 2 == 0:
                    pos = 0
                else:
                    pos = 1
                mtst2.add_state(mts_libs[pos].get_state(mts_name))
            try:
                mt = mtst2.to_motiftree()
                p = mt.to_pose(include_head=1)
            except:
                mtst2.remove_nodes(0)
                continue
            ss = p.secondary_structure()
            seq = p.designable_sequence()
            new_seq = seq[7:-8]
            new_ss = ss[7:-8]
            seq = prime5_seq + new_seq + prime3_seq
            ss  = prime5_ss  + new_ss  + prime3_ss
            fseq = ""
            for e in seq:
                if e == 'T':
                    fseq += 'U'
                else:
                    fseq += e
            seq = fseq
            #print seq
            #print ss
            best = 0
            best_seq = None
            fail = 0
            if len(p.chains()) > 1:
                mtst2.remove_nodes(0)
                continue

            for i in range(1):
                try:
                    results = designer.design(ss, seq)
                    score = results[0]['end'][2]['finalscore']
                    d_new_seq = results[0]['end'][0]
                except:
                    fail = 1
                    print seq
                    print ss
                    p.to_pdb('test.pdb')
                    exit()
                    break
                if score > best:
                    best_seq = d_new_seq
            if not fail:
                print score, best_seq, ss
                f.write(str(count) + " " + str(score) + " " + best_seq + " " + ss + "\n")
            mtst2.remove_nodes(0)
        f.close()


    def get_proper_sequence(self, seq, ss):

        prime5_seq = 'CTAGGATATGG'
        prime3_seq = 'CCTAAGTCCTAG'
        prime5_ss  = '(((((((..(('
        prime3_ss  = '))...)))))))'

        new_seq = seq[7:-8]
        new_ss = ss[7:-8]
        seq = prime5_seq + new_seq + prime3_seq
        ss  = prime5_ss  + new_ss  + prime3_ss
        fseq = ""
        for e in seq:
            if e == 'T':
                fseq += 'U'
            else:
                fseq += e
        seq = fseq
        return seq, ss




    def test_tecto_2(self):
        seq, ss = get_wc_flow()
        flow_ss_tree = secondary_structure_tree.SecondaryStructureTree(ss, seq)
        chip_ss_tree = secondary_structure_tree.SecondaryStructureTree(ss, seq)

        met = get_met(flow_ss_tree, chip_ss_tree)
        mtst = met.get_mtst()
        p = mtst.to_pose()
        chains = p.chains()
        chains.pop()
        basepairs = []
        res = []
        for c in chains:
            res.extend(c.residues)
        for bp in p.basepairs:
            if bp.res1 in res and bp.res2 in res:
                basepairs.append(bp)

        m = motif.Motif()
        m.structure._build_chains(res)
        m.structure._cache_coords()
        m.basepairs = basepairs
        m._cache_basepair_frames()
        m.setup_basepair_ends()

        mts = motif_tree_state.motif_to_state(m)
        mtst2 = motif_tree_state.MotifTreeStateTree(mts)

        twoways = motif_tree_state.MotifTreeStateLibrary(motif_type.TWOWAY)
        helixs = motif_tree_state.MotifTreeStateLibrary(motif_type.HELIX)
        mts_libs = [helixs, twoways]

        f = open('/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/motif_tree_path_refiner/out')
        lines = f.readlines()
        f.close()
        lines.pop(0)

        mtspo = rnamake.motif_tree_state_path_optimizer.MotifTreeStatePathOptimizer()
        designer = sequence_designer.SequenceDesigner()

        f = open('sequences_min.txt')
        seq_lines = f.readlines()
        f.close()

        allowed_counts = []
        for l in seq_lines:
            spl = l.split()
            allowed_counts.append(int(spl[0]))

        f = open('sequences_final.txt', 'w')

        pos = 0
        count = -1
        for l in lines:
            count += 1
            if count not in allowed_counts:
                continue
            spl = l.split('|')
            mts_spl = spl[1].split()
            mts_spl.pop(0)
            for i, mts_name in enumerate(mts_spl):
                if i % 2 == 0:
                    pos = 0
                else:
                    pos = 1
                mtst2.add_state(mts_libs[pos].get_state(mts_name))
            try:
                mt = mtst2.to_motiftree()
                p = mt.to_pose(include_head=1)
            except:
                mtst2.remove_nodes(0)
                continue

            results = mtspo.optimize(mtst2)

            all_mtsts = [mtst2]

            if len(results) > 0:
                all_mtsts.extend(results)


            for mtst in all_mtsts:
                try:
                    mt = mtst2.to_motiftree()
                    p = mt.to_pose(include_head=1)
                except:
                    continue

                ss = p.secondary_structure()
                seq = p.designable_sequence()

                seq, ss = self.get_proper_sequence(seq, ss)

                best = 0
                best_seq = None
                fail = 0
                if len(p.chains()) > 1:
                    mtst2.remove_nodes(0)
                    continue

                for i in range(10):
                    try:
                        results = designer.design(ss, seq)
                        score = results[0]['end'][2]['finalscore']
                        d_new_seq = results[0]['end'][0]
                    except:
                        break
                    if score > best:
                        best_seq = d_new_seq
                if not fail:
                    print count, score, best_seq, ss
                    f.write(str(count) + " " + str(score) + " " + best_seq + " " + ss + "\n")
            mtst2.remove_nodes(0)
        f.close()






def main():
    unittest.main()

if __name__ == '__main__':
    main()
