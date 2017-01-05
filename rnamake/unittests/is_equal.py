import numerical

def are_atom_equal(a1, a2, threshold=0.00001):
    return numerical.are_points_equal(a1.coords, a2.coords, threshold) and a1.name == a2.name


def are_atoms_equal(atoms1, atoms2):
    if len(atoms1) != len(atoms2):
        raise ValueError("cannot call are_atoms_equal with atom lists of different sizes")

    for i in range(len(atoms1)):
        result = are_atom_equal(atoms1[i], atoms2[i])
        if not result:
            return 0

    return 1


def are_residues_equal(r1, r2, check_uuid=1, threshold=0.00001):
    if r1.name != r2.name:
        return 0

    if r1.uuid != r2. uuid and check_uuid:
        return 0

    for i, a in enumerate(r1):
        name = ""
        if a is not None:
            name = a.name
        if a is None and r2.has_atom(index=i):
            return 0
        if r2.has_atom(index=i) == False and a is not None:
            return 0
        if a is None and r2.has_atom(index=i) == False:
            continue

        result = are_atom_equal(a, r2.get_atom(name), threshold)
        if not result:
            return 0

    return 1


def are_chains_equal(c1, c2, check_uuid=1):

    if len(c1) != len(c2):
        return 0

    for i in range(len(c1)):
        result = are_residues_equal(c1.residue(i), c2.residue(i), check_uuid)
        if not result:
            return 0
    return 1


def are_structure_equal(s1, s2, check_uuid=1):
    if len(s1) != len(s2):
        return 0

    for i in range(len(s1)):
        result = are_chains_equal(s1.chain(i), s2.chain(i), check_uuid)
        if not result:
            return 0

    return 1
