import random
from rnamake import secondary_structure

def fill_basepairs_in_ss(ss):
    pairs = ["AU", "UA", "GC", "CG"]
    for bp in ss.basepairs:
        if bp.res1.name == "N" and bp.res2.name == "N":
            p = random.choice(pairs)
            bp.res1.name = p[0]
            bp.res2.name = p[1]

    for m in ss.motifs:
        for i, end in enumerate(m.ends):
            m.end_ids[i] = secondary_structure.assign_end_id_new(m, end)