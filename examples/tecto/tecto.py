import rnamake
import random


def get_first_tecto(attempt=1):
    """
    constructs a RNA that includes two tetraloop receptor tertiary contacts
    seperated between 5 to 10 random helical and twoway junction motifs
    """
    clogger = rnamake.logger.get_logger("tecto:get_first_tecto()")
    clogger.info("Attempting to setup first tetco, try: " + str(attempt))
    gaaa_tetraloop = rnamake.motif.Motif("GAAA_tetraloop")
    bp1 = gaaa_tetraloop.ends[0].uuid
    bp2 = gaaa_tetraloop.ends[1].uuid
    mlibs = [ rnamake.motif_library.ideal_helix_lib(),
              rnamake.motif_library.unique_twoway_lib() ]
    mt = rnamake.motif_tree.MotifTree()
    mt.add_motif(gaaa_tetraloop, end_index=2)
    size = random.randrange(5,10)
    count = 0
    pos = 0
    i = 0
    while i < size:
        if i % 2 == 0:
            pos = 0
        else:
            pos = 1
        m = random.choice(mlibs[pos].motifs())
        node = mt.add_motif(m)
        if node:
            i += 1
        count += 1
        if count > 1000:
            break

    node = mt.add_motif(gaaa_tetraloop, end_index=0)
    #could not add last tetraloop receptor something went wrong
    #recall function
    if node is None:
        return get_first_tecto(attempt+1)
    return mt.to_pose(chain_closure=1), bp1, bp2


def get_second_tecto(self, p, uuid1, uuid2):
    pass

p, uuid1, uuid2 = get_first_tecto()
p.to_pdb()
exit()
bp1s = p.get_basepair(bp_uuid=uuid1)
bp2s = p.get_basepair(bp_uuid=uuid2)
bp1 = bp1s[0]
bp2 = bp2s[1]
chains = []
basepairs = []
for c in p.chains():
    for bp in (bp1,bp2):
        for r in bp.residues():
            if c.first() == r or c.last() == r:
                if c not in chains:
                    chains.append(c)
res = []
for c in chains:
    p.structure.chains.remove(c)
    res.extend(c.residues)
for bp in p.basepairs:
    if bp.res1 in res and bp.res2 in res:
        basepairs.append(bp)
for bp in basepairs:
    p.basepairs.remove(bp)
p.structure._cache_coords()
p._cache_basepair_frames()

m = rnamake.motif.Motif()
m.structure.chains.extend(chains)
m.structure._cache_coords()
m.basepairs = basepairs

m.to_pdb("start2.pdb")
p.to_pdb("start.pdb")
beads = m.get_beads(m.ends)

sl = rnamake.steric_lookup.StericLookup(p)
sl.add_beads(beads)
mtss = rnamake.build_task_astar.MotifTreeStateSearch(accept_score=15,
                                                     max_solutions=1,
                                                     verbose=1)
solutions = mtss.search(bp1.state(), bp2.state(), lookup=sl)
mtst = solutions[0].to_mtst()
mtst.to_pdb("solution.pdb")








