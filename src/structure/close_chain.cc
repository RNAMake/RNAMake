//
// Created by Joseph Yesselman on 4/4/18.
//

#include <math.h>
#include "structure/close_chain.h"

namespace structure {

math::Matrix
create_coord_system(
        AtomOPs const & atoms) {
    auto e1 = (atoms[0]->coords() - atoms[1]->coords()).normalize();
    auto e3 = cross(e1, atoms[2]->coords() - atoms[1]->coords()).normalize();
    auto e2 = cross(e3, e1);
    auto m = math::Matrix(
            e1.x(), e1.y(), e1.z(),
            e2.x(), e2.y(), e2.z(),
            e3.x(), e3.y(), e3.z());
    m = m.transpose();
    return m;
}

math::Vector
m_dot_v(
        math::Matrix const & m,
        math::Vector const & v) {
    return math::Vector(
            m.xx() * v.x() + m.xy() * v.y() + m.xz() * v.z(),
            m.yx() * v.x() + m.yy() * v.y() + m.yz() * v.z(),
            m.zx() * v.x() + m.zy() * v.y() + m.zz() * v.z());
}

AtomOP
virtual_atom(
        String const & name,
        float l,
        float theta,
        float phi,
        AtomOPs const & parent_atoms) {
    theta = to_radians(theta);
    phi = to_radians(phi);
    auto m = create_coord_system(parent_atoms);
    auto v = math::Vector(
            l * cos(theta),
            l * sin(theta) * cos(phi),
            l * sin(theta) * sin(phi));
    auto coords = m_dot_v(m, v);
    coords += parent_atoms[0]->coords();
    return std::make_shared<Atom>(name, coords);
}

float
to_radians(
        float degrees) {
    return (degrees * M_PI) / 180;
}

math::Vector
get_projection(
        math::Point const & coord,
        math::Point const & current_pos,
        math::Vector const & projection_axis) {
    auto r = coord - current_pos;
    return r - r.dot(projection_axis) * projection_axis;
}

math::Matrix
axis_angle_to_rot_matrix(
        float angle,
        math::Vector const & al) {

    auto c = cos(angle);
    auto s = sin(angle);
    auto t = 1.0 - c;
    auto al_norm = al;
    al_norm = al_norm.normalize();
    auto m00 = c + al_norm.x() * al_norm.x() * t;
    auto m11 = c + al_norm.y() * al_norm.y() * t;
    auto m22 = c + al_norm.z() * al_norm.z() * t;

    auto tmp1 = al_norm.x() * al_norm.y() * t;
    auto tmp2 = al_norm.z() * s;
    auto m10 = tmp1 + tmp2;
    auto m01 = tmp1 - tmp2;
    tmp1 = al_norm.x() * al_norm.z() * t;
    tmp2 = al_norm.y() * s;
    auto m20 = tmp1 - tmp2;
    auto m02 = tmp1 + tmp2;
    tmp1 = al_norm.y() * al_norm.z() * t;
    tmp2 = al_norm.x() * s;
    auto m21 = tmp1 + tmp2;
    auto m12 = tmp1 - tmp2;

    return math::Matrix(m00, m01, m02,
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
    auto x = math::Vector(matrix.xx(), matrix.xy(), matrix.xz());
    auto current_atom_xyz = parent_atoms[0]->coords();
    auto weighted_sine = 0.0;
    auto weighted_cosine = 0.0;
    for (int i = 0; i < match_atoms_1.size(); i++) {
        auto rho1 = get_projection(match_atoms_1[i]->coords(), current_atom_xyz, x);
        auto rho2 = get_projection(match_atoms_2[i]->coords(), current_atom_xyz, x);
        auto current_sine = which_dir * dot_product(x, cross(rho1, rho2));
        auto current_cosine = dot_product(rho1, rho2);
        weighted_sine += current_sine;
        weighted_cosine += current_cosine;
    }

    auto twist_torsion = atan2(weighted_sine, weighted_cosine);
    auto m = axis_angle_to_rot_matrix(twist_torsion, x);
    for (auto & daughter_atom : daughter_atoms) {
        auto new_coords = m_dot_v(m, daughter_atom->coords() - current_atom_xyz) + current_atom_xyz;
        daughter_atom->coords(new_coords);
    }

}

math::Matrix
get_res_ref_frame(
        ResidueOP r) {
    auto vec1 = math::Point();
    auto vec2 = math::Point();
    auto beads = r->get_beads();
    auto b = beads[0];
    for (auto const & bead : beads) {
        if (bead.btype() == structure::BeadType::BASE) {
            b = bead;
            break;
        }
    }

    if (r->name() == "A" || r->name() == "G") {
        vec1 = (r->get_atom("N9")->coords() - r->get_atom("C1'")->coords()).normalize();
        vec2 = (r->get_atom("N9")->coords() - b.center()).normalize();
    } else {
        vec1 = (r->get_atom("N1")->coords() - r->get_atom("C1'")->coords()).normalize();
        vec2 = (r->get_atom("N1")->coords() - b.center()).normalize();
    }
    auto cross = vec1.cross(vec2);
    auto m = math::Matrix(
            vec1.x(), vec1.y(), vec1.z(),
            vec2.x(), vec2.y(), vec2.z(),
            cross.x(), cross.y(), cross.z());
    m.unitarize();
    return m;

}


void
replace_missing_phosphate_backbone(
        ResidueOP r,
        ResidueOP r_template) {

    auto ref_frame_1 = get_res_ref_frame(r);
    auto ref_frame_2 = get_res_ref_frame(r_template);

    auto rot = math::Matrix();
    dot(ref_frame_1.transpose(), ref_frame_2, rot);
    auto trans = -r_template->center();
    auto c4p_atom = r->get_atom("C4'");

    auto t = math::Transform(rot, trans);

    r_template->transform(t);
    r_template->move(c4p_atom->coords() - r_template->get_atom("C4'")->coords());

    auto atom_names = Strings{"C5'", "O5'", "P", "OP1", "OP2"};
    auto atoms = AtomOPs();
    for (auto const & name : atom_names) {
        auto a = r_template->get_atom(name);
        atoms.push_back(std::make_shared<Atom>(a->name(), a->coords()));
    }

    for (auto const & a : r->atoms()) {
        if (a == nullptr) { continue; }
        if (std::find(atom_names.begin(), atom_names.end(), a->name()) != atom_names.end()) { continue; }
        atoms.push_back(a);
    }

    r->setup_atoms(atoms);

}

void
close_chain(
        ChainOP chain) {
    if(chain == nullptr) {
        return; 
    }
    auto r_template = ResidueOP(nullptr);
    for (auto const & r : chain->residues()) {
        if(r == nullptr)  {
            continue;
        }
        if (r->get_atom("P") != nullptr && r->get_atom("OP1") != nullptr && r->get_atom("OP2") != nullptr) {
            r_template = std::make_shared<Residue>(*r);
            break;
        }
    }

    if (r_template == nullptr) {
        std::cout << " cannot rebuild phosphates, there is no residue with all phosphate atoms" << std::endl;
        return;
    }
    for (int i = 0; i < chain->length() - 1; i++) {
        auto res1 = chain->residues()[i];
        auto res2 = chain->residues()[i + 1];

        if (res1->get_atom("P") == nullptr) {
            replace_missing_phosphate_backbone(res1, r_template);
        }
        if (res2->get_atom("P") == nullptr) {
            replace_missing_phosphate_backbone(res2, r_template);
        }

        if (res1->connected_to(*res2, 2.0)) {
            continue;
        }

        auto res1_atoms = AtomOPs();
        auto res2_atoms = AtomOPs();
        auto res1_atom_names = Strings{"C4'", "C3'", "O3'"};
        auto res2_atom_names = Strings{"C5'", "O5'", "OP2", "OP1"};

        auto fail = 0;
        for (auto const & a_name : res1_atom_names) {
            auto a = res1->get_atom(a_name);
            if (a == nullptr) {
                fail = 1;
                break;
            }
            res1_atoms.push_back(a);
        }
        for (auto const & a_name : res2_atom_names) {
            auto a = res2->get_atom(a_name);
            if (a == nullptr) {
                fail = 1;
                break;
            }
            res2_atoms.push_back(a);
        }
        if (fail) { continue; }

        auto ovl1 = virtual_atom("OVL1", 1.606497, 60.314519, 0.0,
                                 AtomOPs{res1->get_atom("O3'"), res1->get_atom("C3'"), res1->get_atom("C4'")});
        auto ovl2 = virtual_atom("OVL2", 1.593180, 71.059360, 0.0,
                                 AtomOPs{ovl1, res1->get_atom("O3'"), res1->get_atom("C3'")});
        auto ovu1 = virtual_atom("OVU1", 1.593103, 71.027062, 114.600417,
                                 AtomOPs{res2->get_atom("P"), res2->get_atom("O5'"), res2->get_atom("OP2")});


        auto match_atoms_1 = AtomOPs {res1->get_atom("O3'"), ovl1, ovl2};
        auto match_atoms_2 = AtomOPs {ovu1, res2->get_atom("P"), res2->get_atom("O5'")};

        for (int j = 0; j < 100; j++) {
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
                    AtomOPs{ovu1, res2->get_atom("OP1"), res2->get_atom("OP2"), res2->get_atom("P"),
                            res2->get_atom("O5'")},
                    match_atoms_1, match_atoms_2);

        }


    }

}

}
