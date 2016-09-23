import rnamake.motif_factory as mf
import rnamake.resource_manager as rm
import rnamake.motif_type as motif_type
import rnamake.motif as motif

from rnamake import util

import os

file_path = util.base_dir(os.path.realpath(__file__))

##############################################################################
# LOADING MOTIFS                                                             #
##############################################################################
# motifs are the main objects of rnamake. They can be loaded in multiple ways
# from a preformated motif directory
path = file_path

m1 = mf.factory.motif_from_file(path+"resources/TWOWAY.1GID.4")

#from a pdb
m2 = mf.factory.motif_from_file(path+"resources/motif.pdb")

#from the resource manager
m3 = rm.manager.get_motif(name="HELIX.IDEAL")

#the resource manager can load any motif that is the libraries that rnamake
#contains. It can also get other resources that we will dicuss later.

##############################################################################
# BASIC MOTIF FEATURES                                                       #
##############################################################################
# every object in rnamake that contains structural information can printed
# out in pdb format, make sure to include to .pdb at the end
m1.to_pdb("test_output.pdb")

# In addition to pdb output, all objects including motif can be formatted to
# a text string. This allows for storage of motifs in text files and also
# in sqlite databases
s = m1.to_str()
mnew = motif.str_to_motif(s)

#take a look at the format if your are curious
f = open("motif.txt", "w")
f.write(s)
f.close()

#possibly the most important feature of motifs is that they track the
#basepair ends of each motif in .ends variable. These ends are how motifs
#are combined together by overlapping these edges

print "m.ends => basepair ends of motifs"
for bp in m1.ends:
    print bp.name()

print "##############################################################################"
print "m.chains() => gets all residues in motif"
for c in m1.chains():
    print c


print "##############################################################################"
print "m.residues() => gets all residues in motif"
for r in m1.residues():
    print r
print "##############################################################################"


print "m.sequence() => ", m1.sequence()
print "m.secondary_structure() => ", m1.dot_bracket()
print "m.mtype => ", motif_type.type_to_str(m3.mtype)

m1 = rm.manager.get_motif(name="TWOWAY.1GID.4")

print "##############################################################################"
print "aligning one motif to another"

motif.align_motif(m3.ends[1].state(), m1.ends[0], m1)
m1.to_pdb("aligned.pdb")
m3.to_pdb("aligned_to.pdb")


print "##############################################################################"
print "expected files: "
print " test_output.pdb =>  example of printing out of a motif"
print " aligned.pdb => the motif aligned by its end"
print " aligned_to.pdb => motif that served as reference for aligned.pdb to align to"


