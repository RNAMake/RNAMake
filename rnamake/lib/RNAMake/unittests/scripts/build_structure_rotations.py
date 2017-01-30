from rnamake import structure, settings, residue_type, basic_io
from rnamake.unittests.instances import transform_instances

path = settings.RESOURCES_PATH + "/start/start.pdb"
out_path = settings.LIB_PATH + "/lib/RNAMake/unittests/unittest_resources/structure/"

f = open(out_path+"structure_transformations.dat", "w")
rts = residue_type.ResidueTypeSet()

for i in range(10):
    s = structure.structure_from_pdb(path, rts)
    t = transform_instances.transform_random()

    s.transform(t)
    f.write(basic_io.matrix_to_str(t.rotation()) + "\n")
    f.write(basic_io.point_to_str(t.translation()) + "\n")
    f.write(s.to_str() + "\n")

f.close()

