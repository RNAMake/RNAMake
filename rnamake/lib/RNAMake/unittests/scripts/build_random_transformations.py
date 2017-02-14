from rnamake import residue, settings, residue_type, basic_io
from rnamake.unittests.instances import transform_instances

path = settings.UNITTEST_PATH + "resources/res_strs.dat"
f = open(path)
lines = f.readlines()
f.close()

out_path = settings.LIB_PATH + "/lib/RNAMake/unittests/unittest_resources/math/"

f = open(out_path+"random_transformations.dat", "w")

for i in range(100):
    t = transform_instances.transform_random()

    f.write(basic_io.matrix_to_str(t.rotation()) + "|")
    f.write(basic_io.point_to_str(t.translation()) + "\n")
f.close()

