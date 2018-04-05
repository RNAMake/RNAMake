//
// Created by Joseph Yesselman on 4/4/18.
//

#include <math.h>
#include "structure/close_chain.h"

Matrix
create_coord_system(
        AtomOPs const & atoms) {
    auto e1 = (atoms[0]->coords() - atoms[1]->coords()).normalize();
    auto e3 = cross(e1, atoms[2]->coords() - atoms[1]->coords()).normalize();
    auto e2 = cross(e3, e1);
    auto m = Matrix(
                 e1.x(), e1.y(), e1.z(),
                 e2.x(), e2.y(), e2.z(),
                 e3.x(), e3.y(), e3.z());
    m = m.transpose();
    return m;
}

Vector
m_dot_v(
        Matrix const & m,
        Vector const & v) {
    return Vector(
            m.xx()*v.x() + m.xy()*v.y() + m.xz()*v.z(),
            m.yx()*v.x() + m.yy()*v.y() + m.yz()*v.z(),
            m.zx()*v.x() + m.zy()*v.y() + m.zz()*v.z());
}

AtomOP
virtual_atom(
        String const & name,
        float l,
        float theta,
        float phi,
        AtomOPs const & parent_atoms) {
    theta = to_radians(theta);
    phi   = to_radians(phi);
    auto m = create_coord_system(parent_atoms);
    auto v = Vector(
            l*cos(theta),
            l*sin(theta)*cos(phi),
            l*sin(theta)*sin(phi));
    auto coords = m_dot_v(m, v);
    coords += parent_atoms[0]->coords();
    return std::make_shared<Atom>(name, coords);
}

float
to_radians(
        float degrees) {
    return (degrees*M_PI)/180;
}

Vector
get_projection(
        Point const & coord,
        Point const & current_pos,
        Vector const & projection_axis) {
    auto r = coord - current_pos;
    return r - r.dot(projection_axis)*projection_axis;
}

Matrix
axis_angle_to_rot_matrix(
        float angle,
        Vector const & al) {

    auto c = cos(angle);
    auto s = sin(angle);
    auto t = 1.0 - c;
    auto al_norm = al;
    al_norm = al_norm.normalize();
    auto m00 = c + al_norm.x()*al_norm.x()*t;
    auto m11 = c + al_norm.y()*al_norm.y()*t;
    auto m22 = c + al_norm.z()*al_norm.z()*t;

    auto tmp1 = al_norm.x()*al_norm.y()*t;
    auto tmp2 = al_norm.z()*s;
    auto m10 = tmp1 + tmp2;
    auto m01 = tmp1 - tmp2;
    tmp1 = al_norm.x()*al_norm.z()*t;
    tmp2 = al_norm.y()*s;
    auto m20 = tmp1 - tmp2;
    auto m02 = tmp1 + tmp2;
    tmp1 = al_norm.y()*al_norm.z()*t;
    tmp2 = al_norm.x()*s;
    auto m21 = tmp1 + tmp2;
    auto m12 = tmp1 - tmp2;

    return Matrix(m00, m01, m02,
                  m10, m11, m12,
                  m20, m21, m22);

}

void
close_torsion(
        int which_dir,
        AtomOPs const & parent_atoms,
        AtomOPs daughter_atoms,
        AtomOPs const & match_atoms_1,
        AtomOPs const & match_atoms_2) {
    auto matrix = create_coord_system(parent_atoms);
    matrix = matrix.transpose();
    auto x = Vector(matrix.xx(), matrix.xy(), matrix.xz());
    auto current_atom_xyz = parent_atoms[0]->coords();
    auto weighted_sine = 0.0;
    auto weighted_cosine = 0.0;
    for(int i = 0; i < match_atoms_1.size(); i++) {
        auto rho1 = get_projection(match_atoms_1[i]->coords(), current_atom_xyz, x);
        auto rho2 = get_projection(match_atoms_2[i]->coords(), current_atom_xyz, x);
        auto current_sine = which_dir * dot_product(x, cross(rho1, rho2));
        auto current_cosine = dot_product(rho1, rho2);
        weighted_sine += current_sine;
        weighted_cosine += current_cosine;
    }

    auto twist_torsion = atan2(weighted_sine, weighted_cosine);
    auto m = axis_angle_to_rot_matrix(twist_torsion, x);
    for(auto & daughter_atom : daughter_atoms) {
        auto new_coords = m_dot_v(m, daughter_atom->coords() - current_atom_xyz) + current_atom_xyz;
        daughter_atom->coords(new_coords);
    }

}

void
close_chain(
        ChainOP chain) {

    for(int i = 0; i < chain->length()-1; i++) {

        auto res1 = chain->residues()[i];
        auto res2 = chain->residues()[i+1];

        //if(res1->connected_to(*res2, 2.0)) {
        //    continue;
        //}

        auto res1_atoms = AtomOPs();
        auto res2_atoms = AtomOPs();
        auto res1_atom_names = Strings{"C4'", "C3'", "O3'"};
        auto res2_atom_names = Strings{"C5'","O5'","OP2","OP1"};

        auto fail = 0;
        for(auto const & a_name : res1_atom_names) {
            auto a = res1->get_atom(a_name);
            if(a == nullptr) { fail = 1; break; }
            res1_atoms.push_back(a);
        }
        for(auto const & a_name : res2_atom_names) {
            auto a = res2->get_atom(a_name);
            if(a == nullptr) { fail = 1; break; }
            res2_atoms.push_back(a);
        }
        if(fail) { continue; }

        auto ovl1 = virtual_atom("OVL1", 1.606497, 60.314519, 0.0,
                                 AtomOPs{res1->get_atom("O3'"), res1->get_atom("C3'"), res1->get_atom("C4'")});
        auto ovl2 = virtual_atom("OVL2", 1.593180, 71.059360, 0.0,
                                 AtomOPs{ovl1, res1->get_atom("O3'"), res1->get_atom("C3'")});
        auto ovu1 = virtual_atom("OVU1", 1.593103, 71.027062, 114.600417,
                                 AtomOPs{res2->get_atom("P"), res2->get_atom("O5'"), res2->get_atom("OP2")});


        auto match_atoms_1 = AtomOPs { res1->get_atom("O3'"), ovl1,                ovl2                  };
        auto match_atoms_2 = AtomOPs { ovu1,                  res2->get_atom("P"), res2->get_atom("O5'") };

        for(int j = 0; j < 100; j++) {
            close_torsion(
                    +1,
                    AtomOPs{res1->get_atom("O3'"), res1->get_atom("C3'"), res1->get_atom("C4'")},
                    AtomOPs{ovl1, ovl2},
                    match_atoms_1, match_atoms_2);

            close_torsion(
                    +1,
                    AtomOPs{ovl1, res1->get_atom("O3'"), res1->get_atom("C3'")},
                    AtomOPs{ovl2},
                    match_atoms_1, match_atoms_2);

            close_torsion(
                    -1,
                    AtomOPs{res2->get_atom("P"), res2->get_atom("O5'"), res2->get_atom("C5'")},
                    AtomOPs{ovu1, res2->get_atom("OP1"), res2->get_atom("OP2")},
                    match_atoms_1, match_atoms_2);

            close_torsion(
                    -1,
                    AtomOPs{res2->get_atom("O5'"), res2->get_atom("C5'"), res2->get_atom("C4'")},
                    AtomOPs{ovu1, res2->get_atom("OP1"), res2->get_atom("OP2"), res2->get_atom("P")},
                    match_atoms_1, match_atoms_2);

            close_torsion(
                    -1,
                    AtomOPs{res2->get_atom("C5'"), res2->get_atom("C4'"), res2->get_atom("C3'")},
                    AtomOPs{ovu1, res2->get_atom("OP1"), res2->get_atom("OP2"), res2->get_atom("P"), res2->get_atom("O5'")},
                    match_atoms_1, match_atoms_2);

        }


    }

}