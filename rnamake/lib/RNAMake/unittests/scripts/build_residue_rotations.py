from rnamake import residue, settings, residue_type, basic_io
from rnamake.unittests.instances import transform_instances

path = settings.UNITTEST_PATH + "resources/res_strs.dat"
f = open(path)
lines = f.readlines()
f.close()

out_path = settings.LIB_PATH + "/lib/RNAMake/unittests/unittest_resources/residue/"

f = open(out_path+"residue_transformations.dat", "w")
rts = residue_type.ResidueTypeSet()

for i in range(10):
    r = residue.Residue.from_str(lines[0], rts)
    t = transform_instances.transform_random()

    r.transform(t)
    print t.rotation()
    f.write(basic_io.matrix_to_str(t.rotation()) + "\n")
    f.write(basic_io.point_to_str(t.translation()) + "\n")
    f.write(r.to_str() + "\n")
    break
f.close()

