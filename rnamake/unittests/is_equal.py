import numerical

def are_atom_equal(a1, a2, threshold=0.00001):
    return numerical.are_points_equal(a1.get_coords(), a2.get_coords(), threshold) and \
           a1.get_name() == a2.get_name()


def are_atoms_equal(atoms1, atoms2):
    if len(atoms1) != len(atoms2):
        raise ValueError("cannot call are_atoms_equal with atom lists of different sizes")

    for i in range(len(atoms1)):
        result = are_atom_equal(atoms1[i], atoms2[i])
        if not result:
            return 0

    return 1


def are_residues_equal(r1, r2, check_uuid=1, threshold=0.00001):
    if r1.get_name() != r2.get_name():
        return 0
    if r1.get_chain_id() != r2.get_chain_id():
        return 0
    if r1.get_i_code() != r2.get_i_code():
        return 0
    if r1.get_num() != r2.get_num():
        return 0

    if r1.get_uuid() != r2.get_uuid() and check_uuid:
        return 0

    for i, a in enumerate(r1):
        name = ""
        if a is not None:
            name = a.get_name()
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
        result = are_residues_equal(c1.get_residue(i), c2.get_residue(i), check_uuid)
        if not result:
            return 0
    return 1


def are_structure_equal(s1, s2, check_uuid=1):
    if len(s1) != len(s2):
        return 0

    for i in range(len(s1)):
        result = are_chains_equal(s1.get_chain(i), s2.get_chain(i), check_uuid)
        if not result:
            return 0

    return 1


def are_basepairs_equal(bp1, bp2, check_uuid=1):
    if not numerical.are_points_equal(bp1.get_d(), bp2.get_d()):
        return 0
    if not numerical.are_matrices_equal(bp1.get_r(), bp2.get_r()):
        return 0
    if not numerical.are_points_equal(bp1.get_res1_sugar(), bp2.get_res1_sugar()):
        return 0
    if not numerical.are_points_equal(bp1.get_res2_sugar(), bp2.get_res2_sugar()):
        return 0

    if bp1.get_name() != bp2.get_name():
        return 0
    if bp1.get_bp_type() != bp2.get_bp_type():
        return 0
    if bp1.get_x3dna_bp_type() != bp2.get_x3dna_bp_type():
        return 0

    if check_uuid:
        if bp1.get_uuid() != bp2.get_uuid():
            return 0
        if bp1.get_res1_uuid() != bp2.get_res1_uuid():
            return 0
        if bp1.get_res2_uuid() != bp2.get_res2_uuid():
            return 0

    return 1


def are_rna_strucs_equal(rs1, rs2, check_uuid=1):

    if rs1.get_num_res() != rs2.get_num_res():
        return 0

    for i in range(rs1.get_num_res()):
        result = are_residues_equal(rs1.get_residue(index=i),
                                    rs2.get_residue(index=i),
                                    check_uuid)
        if not result:
            return 0

    if rs1.get_num_basepairs() != rs2.get_num_basepairs():
        return 0

    for i in range(rs1.get_num_basepairs()):
        result = are_basepairs_equal(rs1.get_basepair(index=i),
                                     rs2.get_basepair(index=i),
                                     check_uuid)
        if not result:
            return 0

    if rs1.get_num_ends() != rs2.get_num_ends():
        return 0

    for i in range(rs1.get_num_ends()):
        result = are_basepairs_equal(rs1.get_end(i), rs2.get_end(i), check_uuid)
        if not result:
            return 0

    for i in range(rs1.get_num_ends()):
        if rs1.get_end_id(i) != rs2.get_end_id(i):
            return 0

    return 1