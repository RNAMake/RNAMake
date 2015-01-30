import rnamake
import numpy as np

def test_str_to_motif():
    path = rnamake.settings.UNITTEST_PATH + "/resources/motifs/p4p6"
    m = rnamake.motif.Motif(path)
    f = open("test_str_to_motif.dat", "w")
    f.write( m.to_str() )
    f.close()



test_str_to_motif()
