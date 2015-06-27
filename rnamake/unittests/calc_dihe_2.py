import rnamake
import rnamake.sqlite_library as sqlite_library
import rnamake.motif_factory as mf
import rnamake.util as util
import rnamake.structure
import math
import glob
import numpy as np
import fast_plotting.line_plot as line_plot

def arc_cos(num):
    return math.atan2(math.sqrt(1-num*num), num)

def diff(p1, p2):
    return np.array([p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2]])

def normalize(v):
    mag = util.distance(v, [0,0,0])
    v_norm = [ x / mag for x in v ]
    return v_norm

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

def get_new_dihe_from_pdb(pdb_path):

    dihedral_lists = [
        "C4' O4' C1' C2'".split(),
        "O4' C1' C2' C3'".split(),
        "C1' C2' C3' C4'".split(),
        "C2' C3' C4' O4'".split(),
        "C3' C4' O4' C1'".split()
    ]


    struct = mf.factory.get_structure(pdb_path)

    r = struct.residues()[0]
    cdihes = []
    for dl in dihedral_lists:
        atoms = []
        for name in dl:
            atoms.append(r.get_atom(name))
        cdihes.append(calc_dihe_3(atoms))

    return cdihes

def get_all_dihes():
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

    return all_dihes

def frac_within_cutoff(cdihe, all_dihes, cutoff):
    count = 0
    total = 0
    for dihe in all_dihes:
        diff = 0
        if dihe[2] != 'G':
            continue
        total += 1
        for i in range(4,9):
            diff += (cdihes[i-4] - dihe[i]) ** 2

        diff = math.sqrt(diff / 5)
        if diff < cutoff:
            print dihe
            count += 1

    return float(count) / float(total)



pdb_path = '/Users/josephyesselman/Downloads/3bo4_g.pdb'
cdihes = get_new_dihe_from_pdb(pdb_path)
all_dihes = get_all_dihes()

x = []
y = []

print cdihes
f = open("calc_dihe_data/data_3bo4.dat", "w")

for i in range(0, 100, 1):
    frac = frac_within_cutoff(cdihes, all_dihes, 1)
    exit()
    x.append(i)
    y.append(frac*100)
    f.write(str(i) + " " + str(frac*100) + "\n")
f.close()

fplp = line_plot.FPLinePlot(1, x, y, xlabel="Dihedral RMSD Cutoff ", ylabel="% G Residues",title="Open Complex")
fplp.plot()
