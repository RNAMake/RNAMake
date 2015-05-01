import rnamake
import numpy as np

def test_str_to_chain():
    path = rnamake.settings.UNITTEST_PATH + "/resources/motifs/p4p6"
    m = rnamake.motif.Motif(path)
    f = open("test_str_to_chain.dat", "w")
    for c in m.chains():
        f.write(c.to_str() + "\n")
    f.close()

def test_str_to_structure():
    path = rnamake.settings.UNITTEST_PATH + "/resources/motifs/p4p6"
    m = rnamake.motif.Motif(path)
    f = open("test_str_to_structure.dat", "w")
    f.write( m.structure.to_str() )
    f.close()

def test_transform():
    path = rnamake.settings.UNITTEST_PATH + "/resources/motifs/p4p6"
    m = rnamake.motif.Motif(path)
    r = np.random.random([3,3])
    d = np.random.random([3])
    t = rnamake.transform.Transform(r, d)
    m.structure.transform(t)
    f = open("test_transform.dat", "w")
    f.write(rnamake.basic_io.matrix_to_str(r) + "\n")
    f.write(rnamake.basic_io.point_to_str(d) + "\n")
    f.write(m.structure.to_str())
    m.structure.to_pdb("transformed3.pdb")
    f.close()



#test_str_to_chain()
#test_str_to_structure()
test_transform()
