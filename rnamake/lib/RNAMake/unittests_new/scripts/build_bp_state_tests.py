from rnamake.unittests import instances
from rnamake import settings, basic_io

RES_PATH = settings.LIB_PATH + "/lib/RNAMake/unittests_new/resources/structure/"

def test_get_transforming_r_and_t():

    path = RES_PATH + "get_transforming_r_and_t_test.dat"
    f = open(path, "w")
    for i in range(100):
        bp_state_1 = instances.basepairstate_random()
        bp_state_2 = instances.basepairstate_random()

        r, t = bp_state_1.get_transforming_r_and_t_w_state(bp_state_2)
        f.write(bp_state_1.to_str() + "|" + bp_state_2.to_str() + "|")
        f.write(basic_io.point_to_str(t) + "|")
        f.write(basic_io.matrix_to_str(r) + "\n")

    f.close()

def test_get_transformed_state():
    path = RES_PATH + "test_get_transformed_state.dat"
    f = open(path, "w")
    for i in range(100):
        bp_state_1 = instances.basepairstate_random()
        bp_state_2 = instances.basepairstate_random()

        f.write(bp_state_1.to_str() + "|" + bp_state_2.to_str() + "|")
        r, t = bp_state_1.get_transforming_r_and_t_w_state(bp_state_2)
        t += bp_state_1.d

        new_r, new_d, new_sug = bp_state_2.get_transformed_state(r, t)
        bp_state_2.set(new_r, new_d, new_sug)
        f.write(bp_state_2.to_str() + "\n")
    f.close()

test_get_transforming_r_and_t()
test_get_transformed_state()


