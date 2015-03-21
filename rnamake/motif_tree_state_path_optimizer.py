import itertools
import settings
import motif_tree_state

class MotifTreeStatePathOptimizer(object):
    def __init__(self):
        path = settings.RESOURCES_PATH + "sim_list"
        f = open(path)
        lines = f.readlines()
        f.close()

        lines.pop(0)
        path = settings.PRECOMPUTED_PATH + "/motif_tree_states/TWOWAY_all.new.me"
        self.mts_lib = motif_tree_state.MotifTreeStateLibrary(libpath=path)
        self.sim_dict = {}


        for l in lines:
            spl = l.split()
            key = spl[1]
            sims = spl[4:]
            self.sim_dict[key] = sims

    def optimize(self, mtst):
        mts_alts = []
        multi = 0
        for i, n in enumerate(mtst.nodes[1:]):
            name_spl = n.mts.name.split('-')
            mname = name_spl.pop(0)
            rest = '-'.join(name_spl)
            if mname in self.sim_dict:
                alt_motifs = self.sim_dict[mname]
                if len(alt_motifs) > 1:
                    multi = 1
                all_mts = []
                for alt_name in alt_motifs:
                    mts_name = alt_name + "-" + rest
                    mts = self.mts_lib.get_state(mts_name)
                    all_mts.append(mts)
                mts_alts.append(all_mts)
            else:
                mts_alts.append([n.mts])

        if multi == 0:
            return []


        combos = itertools.product(*mts_alts)
        seq = mtst.to_pose().designable_sequence()


        unique = 0
        mtsts = []
        seen = [seq]
        for i, c in enumerate(combos):
            mtst2 = motif_tree_state.MotifTreeStateTree(mtst.nodes[0].mts)
            fail = 0
            for mts in c:
                node = mtst2.add_state(mts)
                if node is None:
                    fail = 1
                    break
            if fail:
                continue

            seq2 = mtst2.to_pose().designable_sequence()

            if seq2 in seen:
                continue
            seen.append(seq2)

            unique = 1
            mtsts.append(mtst2)

        return mtsts








