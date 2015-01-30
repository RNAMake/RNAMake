import rnamake
import numpy as np
import random

def test_str_to_basepairstate():
    path = rnamake.settings.UNITTEST_PATH + "/resources/motifs/p4p6"
    m = rnamake.motif.Motif(path)
    f = open ("test_str_to_basepairstate.dat", "w")
    for bp in m.basepairs:
        f.write ( bp.state().to_str() + "\n")
    f.close()

def build_get_transforming_r_and_t_test():
    path = rnamake.settings.UNITTEST_PATH + "/resources/motifs/p4p6"
    m = rnamake.motif.Motif(path)
    f = open("get_transforming_r_and_t_test.dat","w")
    for i in range(100):
        bp1 = random.choice(m.basepairs).state()
        bp2 = random.choice(m.basepairs).state()
        r,t = bp1.get_transforming_r_and_t_w_state(bp2)
        f.write(bp1.to_str() + "|" + bp2.to_str() + "|" + \
                rnamake.basic_io.point_to_str(t) + "|" + \
                rnamake.basic_io.matrix_to_str(r) + "\n")
    f.close()



#test_str_to_chain()
test_str_to_basepairstate()
build_get_transforming_r_and_t_test()
