import ensemble_utils
import sys
import re

sequence  = "ggNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNaggauauNNNNNNNNNNagaaggNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNaaagaaacaacaacaacaac"
secstruct = ".....(((((((((((((((....)))))))((((......((((....)))).....))))(((((((....)))))))))))))))...................."

print "material = rna1999"
print "structure s = %s" % (secstruct)

for ii in range(0,len(sequence)):
	print "domain a%d = %s1" % (ii, sequence[ii].upper())

sstr = "s.seq = "
for ii in range(0,len(sequence)):
	sstr += "a%d " % ii

print sstr
print "s.stop = 1.0"
print "prevent = GGGG, CCCC"
print "trials = 10"