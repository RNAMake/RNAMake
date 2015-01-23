import rnamake
import rnamake.unittests.build_task_astar_unittest as build_task_astar_unittest

# In this example we are going to start to use the MotifTreeState search
# which allows us to find combinations of motifs that will get us from
# point A to point B

# Here I have generated a random helix/twoway junction rna
random_tree = build_task_astar_unittest.get_twoway_helix_mts_tree(size=10)
random_tree.to_pdb("prob1.pdb")

# I then grab the first and last node of the tree to define the start
# and end of the search. Both start and end are BasepairState objects. They
# represent the minimal amount of information requred to define a search
# from one basepair orientation to another
start = random_tree.nodes[0].active_states()[0]
end = random_tree.nodes[-1].active_states()[0]
print "starting at:"
print "rotation: ",start.r
print "origin: ",start.d
print "##############################################################################"
print "going to:"
print "rotation: ",end.r
print "origin: ",end.d

mtss = rnamake.build_task_astar.MotifTreeStateSearch(verbose=1,
                                                     max_solutions=1,
                                                     accept_score=20)

print "##############################################################################"
print "building with helices and two way junctions"

solution = mtss.search(start, end)[0]
solution.to_pdb("solution1.pdb")

print "COMPLETE: view problem and solutions, 'pymol prob1.pdb solution1.pdb'"

ns = rnamake.build_task_astar.MotifTreeStateSelector([rnamake.motif_type.TWOWAY],
                                                     "all")

mtss.reset()
mtss.node_selector = ns

print "building with two way junctions"

solution = mtss.search(start, end)[0]
solution.to_pdb("solution2.pdb")

print "COMPLETE: view problem and solutions, 'pymol prob1.pdb solution2.pdb'"

ns = rnamake.build_task_astar.MotifTreeStateSelector([rnamake.motif_type.NWAY])

mtss.reset()
mtss.node_selector = ns

solution = mtss.search(start, end)[0]
solution.to_pdb("solution3.pdb")


