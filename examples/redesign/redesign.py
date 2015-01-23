import rnamake

p = rnamake.pose.Pose("resources/p4p6")
twoway = p.twoways()[3]

segmenter = rnamake.segmenter.Segmenter()
segments = segmenter.apply(p, twoway.ends)
segments.remaining.to_pdb()

start = twoway.ends[0].state()
end = twoway.ends[1].state()
sl = rnamake.steric_lookup.StericLookup()
sl.add_beads(segments.remaining.get_beads(twoway.ends))


mtss = rnamake.build_task_astar.MotifTreeStateSearch(verbose=1,
                                                     max_solutions=1,
                                                     accept_score=7,
                                                     max_node_level=3)

solutions = mtss.search(start, end, lookup=sl)
solutions[0].to_pdb("test.pdb")
