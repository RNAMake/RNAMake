import rnamake.pose_factory as pf
import rnamake.motif_type as motif_type

#A pose is a RNA structure that contains more then one motif. This is not
#rigidly enforced, but does allow for more functionality for larger structures
#like a motif, a pose can be loaded from a preformed directory or from a pdb

p = pf.factory.pose_from_file("resources/p4p6")

print "this structure contains", len(p.motifs(motif_type.TWOWAY)), "twoway junctions"
print "visualize them with pymol twoways.*"
for i, t in enumerate(p.motifs(motif_type.TWOWAY)):
   t.to_pdb("twoways."+str(i)+".pdb")

print "##############################################################################"
print "p.designable_sequence() => get designable sequence"
#print  p.designable_sequence()

print "##############################################################################"
print "p.optimized_sequence() => get optimized sequence"
#print p.optimized_sequence()

for m in p.motifs(motif_type.ALL):
    if m.mtype == motif_type.HELIX:
        continue
    print m.sequence(), m.dot_bracket()




