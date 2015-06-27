import rnamake
import rnamake.sqlite_library as sqlite_library
import rnamake.motif_factory as mf
import rnamake.util as util
import rnamake.structure
import math
import glob
import numpy as np

def arc_cos(num):
    return math.atan2(math.sqrt(1-num*num), num)

def diff(p1, p2):
    return np.array([p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2]])

def normalize(v):
    mag = util.distance(v, [0,0,0])
    v_norm = [ x / mag for x in v ]
    return v_norm

def calc_dihe(atoms):
    d12 =util.distance(atoms[0].coords, atoms[1].coords)
    d23 =util.distance(atoms[1].coords, atoms[2].coords)
    d13 =util.distance(atoms[0].coords, atoms[2].coords)
    d24 =util.distance(atoms[1].coords, atoms[3].coords)
    d34 =util.distance(atoms[2].coords, atoms[3].coords)
    d14 =util.distance(atoms[0].coords, atoms[3].coords)
    p = d12**2 * ( d23**2 + d34**2 - (d24**2)) + \
	            d23**2 * (-(d23**2) + d34**2 + d24**2) + \
	            d13**2 * ( d23**2 - d34**2 + d24**2) - \
	            2 * d23**2 * d14**2

    q = (d12 + d23 + d13) * ( d12 + d23 - d13) * \
		(d12 - d23 + d13) * (-d12 + d23 + d13 ) *\
		(d23 + d34 + d24) * ( d23 + d34 - d24 ) * \
		(d23 - d34 + d24) * (-d23 + d34 + d24 )

    phi = arc_cos(p/math.sqrt(q))
    return phi*180 / math.pi

def calc_dihe_2(atoms):
    ba = diff(atoms[1].coords, atoms[0].coords)
    cb = diff(atoms[2].coords, atoms[1].coords)
    dc = diff(atoms[3].coords, atoms[2].coords)

    cross_ba_cb = np.cross(ba,cb)
    cross_cb_dc = np.cross(cb,dc)
    phi = arc_cos(cross_ba_cb.dot(cross_cb_dc) / (util.distance(cross_ba_cb, [0,0,0])*util.distance(cross_cb_dc, [0,0,0])))
    return phi*180 / math.pi

def calc_dihe_3(atoms):
    b_a = diff(atoms[1].coords, atoms[0].coords)
    b_c = diff(atoms[1].coords, atoms[2].coords)
    c_d = diff(atoms[3].coords, atoms[2].coords)
    b_a = normalize(b_a)
    b_c = normalize(b_c)
    c_d = normalize(c_d)

    n1 = np.cross(b_a, b_c)
    n2 = np.cross(b_c, c_d)
    m  = np.cross(n1, b_c)
    x = np.dot(n1, n2)
    y = np.dot(m, n2)

    return (180 / math.pi) * math.atan2(y, x)


mlib = sqlite_library.MotifSqliteLibrary("ideal_helices")
m = mlib.get('HELIX.IDEAL')
#m.to_pdb('test.pdb')

dihedral_lists = [
    "C4' O4' C1' C2'".split(),
    "O4' C1' C2' C3'".split(),
    "C1' C2' C3' C4'".split(),
    "C2' C3' C4' O4'".split(),
    "C3' C4' O4' C1'".split()
]

path = '/Users/josephyesselman/projects/REDESIGN/resources/non-redundant-rnas/'
files = glob.glob(path+"/*")

prime2 = mf.factory.get_structure('/Users/josephyesselman/Downloads/1u6b_g.pdb')

r = prime2.residues()[0]
cdihes = []
for dl in dihedral_lists:
    atoms = []
    for name in dl:
        atoms.append(r.get_atom(name))
    cdihes.append(calc_dihe_3(atoms))


f = open('all_sugar_dihedrals.txt')
lines = f.readlines()
f.close()

all_dihes = []
for l in lines:
    spl = l.split()
    if len(spl) < 8:
        continue
    spl[4] = float(spl[4])
    spl[5] = float(spl[5])
    spl[6] = float(spl[6])
    spl[7] = float(spl[7])
    spl[8] = float(spl[8])
    all_dihes.append(spl)

cutoff = 25
closest_distance = 100000
closest_hit = 0

worst_distance = 0

for dihe in all_dihes:
    diff = 0
    if dihe[2] != 'G':
        continue
    for i in range(4,9):
        diff += abs(cdihes[i-4] - dihe[i])

    if diff < closest_distance:
        closest_distance = diff
        closest_hit = dihe

    if diff > worst_distance:
        worst_distance = diff

    #if diff < cutoff:
    #    print dihe

print cdihes
#print closest_hit, closest_distance, worst_distance

exit(0)
fsum = open('all_sugar_dihedrals.txt', 'w')

for f in files:
    try:
        print f
        m = rnamake.motif.Motif(f)
    except:
        continue
    print f, len(m.residues())
    res = m.residues()
    for r in res:
        fsum.write(m.name + " " + str(r.num) + " " + str(r.name) + " " + r.chain_id + " ")
        fail = 0
        for dl in dihedral_lists:
            atoms = []
            for name in dl:
                atoms.append(r.get_atom(name))
            try:
                dihe = calc_dihe_3(atoms)
                fsum.write(str(dihe) + " " )
            except:
                break
        fsum.write("\n")


exit(0)
for dl in dihedral_lists:
    atoms = []
    for name in dl:
        atoms.append(m.residues()[0].get_atom(name))
    dihe = calc_dihe(atoms)
    print dihe
