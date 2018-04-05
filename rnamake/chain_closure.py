import math
import numpy as np
import atom
import util
import transform


def create_coord_system(atoms):
    assert(len(atoms) == 3)
    e1 = util.normalize(atoms[0].coords - atoms[1].coords)
    e3 = np.cross(e1, atoms[2].coords - atoms[1].coords)
    e3 = util.normalize(e3)
    e2 = np.cross(e3, e1)
    matrix = np.array([e1, e2, e3]).T
    return matrix


def virtual_atom(name, l, theta, phi, parent_atoms):
    # get radians
    theta = math.radians(theta)
    phi   = math.radians(phi)
    matrix = create_coord_system(parent_atoms)
    v = [l*math.cos(theta),
         l*math.sin(theta)*math.cos(phi),
         l*math.sin(theta)*math.sin(phi)]
    coords = np.dot(matrix, v) + parent_atoms[0].coords
    return atom.Atom(name, coords)


def get_projection(coord, current_pos, projection_axis):
    r = coord - current_pos
    return r - np.dot(r, projection_axis)*projection_axis


def rotation_matrix(angle, direction, point=None):
    """Return matrix to rotate about axis defined by point and direction.

    >>> R = rotation_matrix(math.pi/2, [0, 0, 1], [1, 0, 0])
    >>> np.allclose(numpy.dot(R, [0, 0, 0, 1]), [1, -1, 0, 1])
    True
    >>> angle = (random.random() - 0.5) * (2*math.pi)
    >>> direc = np.random.random(3) - 0.5
    >>> point = np.random.random(3) - 0.5
    >>> R0 = rotation_matrix(angle, direc, point)
    >>> R1 = rotation_matrix(angle-2*math.pi, direc, point)
    >>> is_same_transform(R0, R1)
    True
    >>> R0 = rotation_matrix(angle, direc, point)
    >>> R1 = rotation_matrix(-angle, -direc, point)
    >>> is_same_transform(R0, R1)
    True
    >>> I = np.identity(4, numpy.float64)
    >>> np.allclose(I, rotation_matrix(math.pi*2, direc))
    True
    >>> np.allclose(2, numpy.trace(rotation_matrix(math.pi/2,
    ...                                               direc, point)))
    True

    """
    sina = math.sin(angle)
    cosa = math.cos(angle)
    direction = util.normalize(direction[:3])
    # rotation matrix around unit vector
    R = np.diag([cosa, cosa, cosa])
    R += np.outer(direction, direction) * (1.0 - cosa)
    direction *= sina
    R += np.array([[ 0.0,         -direction[2],  direction[1]],
                      [ direction[2], 0.0,          -direction[0]],
                      [-direction[1], direction[0],  0.0]])
    M = np.identity(4)
    M[:3, :3] = R
    if point is not None:
        # rotation not around origin
        point = np.array(point[:3], dtype=numpy.float64, copy=False)
        M[:3, 3] = point - np.dot(R, point)
    return M


def close_torsion(which_dir, parent_atoms, daughter_atoms, match_atoms_1,
                  match_atoms_2):
    matrix = create_coord_system(parent_atoms)
    x = matrix.T[0]
    current_atom_xyz = parent_atoms[0].coords
    weighted_sine = 0
    weighted_cosine = 0
    for i in range(len(match_atoms_1)):
        rho1 = get_projection(match_atoms_1[i].coords, current_atom_xyz, x)
        rho2 = get_projection(match_atoms_2[i].coords, current_atom_xyz, x)
        current_sine = which_dir * np.dot(x, np.cross(rho1,rho2))
        current_cosine = np.dot(rho1,rho2)
        weighted_sine += current_sine
        weighted_cosine += current_cosine

    twist_torsion = math.atan2(weighted_sine, weighted_cosine)
    for daughter_atom in daughter_atoms:
        daughter_atom.coords = np.dot(rotation_matrix( twist_torsion , x )[:3,:3],
                                      ( daughter_atom.coords - current_atom_xyz )) + \
                               current_atom_xyz

# Copied from Rhiju's Rosetta code: protocols/farna/RNA_LoopCloser.cc
def close_chain(chain):
    r_template = None
    for r in chain.residues:
        if r.get_atom("P") is not None:
            r_template = r
            break

    if r_template is None:
        raise ValueError(
            "cannot close chain no residue has a phostphate")

    r_template = r_template.copy()


    for i in range(0,len(chain.residues)-1):

        res1 = chain.residues[i]
        res2 = chain.residues[i+1]

        if res1.get_atom("P") is None:
            replace_missing_phosphate_backbone(res1, r_template)
        if res2.get_atom("P") is None:
            replace_missing_phosphate_backbone(res2, r_template)


        if res1.connected_to(res2,cutoff=2.0):
            continue

        atoms = []
        res1_atoms = [res1.get_atom(name) for name in [ "C4'", "C3'", "O3'" ]]
        res2_atoms = [res2.get_atom(name) for name in "C5',O5',OP2,OP1".split(",")]
        atoms.extend(res1_atoms)
        atoms.extend(res2_atoms)
        fail=0
        for a in atoms:
            if a is None:
                fail=1
        if fail:
            continue

        ovl1 = virtual_atom("OVL1", 1.606497, 60.314519, 0.0,
                            [res1.get_atom("O3'"), res1.get_atom("C3'"),
                             res1.get_atom("C4'")])

        ovl2 = virtual_atom("OVL2", 1.593180, 71.059360, 0.0,
                            [ovl1, res1.get_atom("O3'"),
                             res1.get_atom("C3'")])
        ovu1 = virtual_atom("OVU1", 1.593103, 71.027062, 114.600417,
                            [res2.get_atom("P"), res2.get_atom("O5'"),
                            res2.get_atom("OP2")])

        res1.atoms.extend([ovl1, ovl2])
        res2.atoms.append(ovu1)

        match_atoms_1 = [res1.get_atom("O3'"), ovl1,               ovl2                ]
        match_atoms_2 = [ovu1,                 res2.get_atom("P"), res2.get_atom("O5'")]

        for j in range(100):
            close_torsion(+1,
                              [res1.get_atom("O3'"), res1.get_atom("C3'"),res1.get_atom("C4'")],
                              [ovl1,ovl2], match_atoms_1, match_atoms_2)

            close_torsion(+1,
                              [ovl1,res1.get_atom("O3'"),res1.get_atom("C3'")],
                              [ovl2], match_atoms_1, match_atoms_2)
            close_torsion(-1,
                              [res2.get_atom("P"),res2.get_atom("O5'"),res2.get_atom("C5'")],
                              [ovu1,res2.get_atom("OP1"),res2.get_atom("OP2")],
                              match_atoms_1, match_atoms_2)
            close_torsion(-1,
                              [res2.get_atom("O5'"),res2.get_atom("C5'"),res2.get_atom("C4'")],
                              [ovu1,res2.get_atom("OP1"),res2.get_atom("OP2"),res2.get_atom("P")],
                              match_atoms_1, match_atoms_2)
            close_torsion(-1,
                              [res2.get_atom("C5'"),res2.get_atom("C4'"),res2.get_atom("C3'")],
                              [ovu1,res2.get_atom("OP1"),res2.get_atom("OP2"),
                               res2.get_atom("P"),res2.get_atom("O5'")],
                              match_atoms_1, match_atoms_2)

        res1.atoms.pop()
        res1.atoms.pop()
        res2.atoms.pop()


def replace_missing_phosphate_backbone(r, r_template):
    ref_frame_1 = get_res_ref_frame(r)
    ref_frame_2 = get_res_ref_frame(r_template)

    rot = util.unitarize(ref_frame_1.T.dot(ref_frame_2))
    t = transform.Transform(rot, -r_template.center())
    r_template.transform(t)
    r_template.move(r.get_atom("C4'").coords - r_template.get_atom("C4'").coords)

    atom_names = ["C5'", "O5'", "P", "OP1", "OP2"]
    new_atoms = []
    for name in atom_names:
        a = r_template.get_atom(name)
        new_atoms.append(atom.Atom(a.name, np.copy(a.coords)))

    for na in new_atoms:
        pos = r.rtype.atom_map[na.name]
        r.atoms[pos] = na


def get_res_ref_frame(r):
    beads =  r.get_beads()
    c = r.center()
    if r.name == "A" or r.name == "G":
        vec1 = util.normalize(np.array(r.get_atom("N9").coords - r.get_atom("C1'").coords))
        vec2 = util.normalize(np.array(r.get_atom("N9").coords - beads[-1].center))
        cross = np.cross(vec1, vec2)
        rot = [vec1, vec2, cross]
    else:
        vec1 = util.normalize(np.array(r.get_atom("N1").coords - r.get_atom("C1'").coords))
        vec2 = util.normalize(np.array(r.get_atom("N1").coords - beads[-1].center))
        cross = np.cross(vec1, vec2)
        rot = [vec1, vec2, cross]

    rot = np.array(rot)
    return util.unitarize(rot)
