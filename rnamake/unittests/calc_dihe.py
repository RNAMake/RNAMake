import rnamake
import rnamake.util as util
import rnamake.structure
import math
import glob

def calc_dihe(atoms):
    d12 =util.distance(atoms[0].coords, atoms[1].coords)
    d23 =util.distance(atoms[1].coords, atoms[2].coords)
    d13 =util.distance(atoms[0].coords, atoms[2].coords)
    d24 =util.distance(atoms[1].coords, atoms[3].coords)
    d34 =util.distance(atoms[2].coords, atoms[3].coords)
    d14 =util.distance(atoms[0].coords, atoms[3].coords)
    p = d12**2 * ( d23**2+ d34**2-d24**2) + \
	            d23**2 * (-d23**2+d34**2+d24**2) + \
	            d13**2 * ( d23**2-d34**2+d24**2) - \
	            2 * d23**2 * d14**2

    q = (d12 + d23 + d13) * ( d12 + d23 - d13) * \
		(d12 - d23 + d13) * (-d12 + d23 + d13 ) *\
		(d23 + d34 + d24) * ( d23 + d34 - d24 ) * \
		(d23 - d34 + d24) * (-d23 + d34 + d24 )

    phi = math.acos(p/math.sqrt(q))
    return phi*180 / math.pi

mlib = rnamake.motif_library.MotifLibrary(rnamake.motif_type.HELIX)
m = mlib.get_motif('HELIX.IDEAL')
#m.to_pdb('test.pdb')

atom_names = "C1' O4' C4' C3'".split()

dihedral_lists = [
    "C4' O4' C1' C2'".split(),
    "O4' C1' C2' C3'".split(),
    "C1' C2' C3' C4'".split(),
    "C2' C3' C4' O4'".split(),
    "C3' C4' O4' C1'".split()
]

path = '/Users/josephyesselman/projects/REDESIGN/resources/non-redundant-rnas/'
files = glob.glob(path+"/*")

prime2 = rnamake.structure.Structure(pdb='/Users/josephyesselman/Downloads/1y0q_g.pdb')

r = prime2.residues()[0]
cdihes = []
for dl in dihedral_lists:
    atoms = []
    for name in dl:
        atoms.append(r.get_atom(name))
    cdihes.append(calc_dihe(atoms))

print cdihes

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

cutoff = 15
for dihe in all_dihes:
    diff = 0
    if dihe[2] != 'G':
        continue
    for i in range(4,9):
        diff += abs(cdihes[i-4] - dihe[i])

    if diff < cutoff:
        print dihe


exit()

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
                dihe = calc_dihe(atoms)
                fsum.write(str(dihe) + " " )
            except:
                break
        fsum.write("\n")

for dl in dihedral_lists:
    atoms = []
    for name in dl:
        atoms.append(m.residues()[0].get_atom(name))
    dihe = calc_dihe(atoms)
    print dihe
