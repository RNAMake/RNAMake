import rnamake

# In this example we are going to redesign a specific part of the P4P6 domain
# structure
p = rnamake.pose.Pose("resources/p4p6")
# get the bounds of a motif we want to redesign
twoway = p.twoways()[3]

#cut out that peices, using the ends as the boundaries to cut
segmenter = rnamake.segmenter.Segmenter()
segments = segmenter.apply(p, twoway.ends)

#visualize the segmentation, this needs to be done so a new peice can be built
#in
segments.remaining.to_pdb("remaining.pdb")
segments.removed.to_pdb("removed.pdb")

#get bounds to setup search
start = twoway.ends[0].state()
end = twoway.ends[1].state()
sl = rnamake.steric_lookup.StericLookup()
sl.add_beads(segments.remaining.get_beads(twoway.ends))


mtss = rnamake.build_task_astar.MotifTreeStateSearch(verbose=1,
                                                     max_solutions=1,
                                                     accept_score=7,
                                                     max_node_level=3)

solutions = mtss.search(start, end, lookup=sl)
s_pose = solutions[0].to_pose()

#put everthing back together
mt = rnamake.motif_tree.MotifTree(segments.remaining, sterics=0)
parent_end = mt.nodes[0].motif.get_basepair(bp_uuid=twoway.ends[0].uuid)[0]
parent_end2 = mt.nodes[0].motif.get_basepair(bp_uuid=twoway.ends[1].uuid)[0]
node = rnamake.motif_tree.MotifTreeNode(s_pose, 1, 1, 0)
rnamake.motif_tree.MotifTreeConnection(mt.nodes[0], node, parent_end, s_pose.ends[1])
mt.nodes.append(node)
new_pose = mt.to_pose(include_head=1,chain_closure=1)
new_pose.to_pdb("new_design.pdb")

print "sequence and secondary structure for new construct, also can visual it as"
print "new_design.pdb"

print new_pose.secondary_structure()
print new_pose.optimized_sequence()






