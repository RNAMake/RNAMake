//
//  residue.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <sstream>
#include <fstream>

//RNAMake Headers
#include <structure/residue.h>

namespace structure {

    math::Point
    center(
            Atoms const &atoms) {
        auto center = math::Point();
        for (auto const &a : atoms) {
            center += a.get_coords();
        }
        return center / float(atoms.size());
    }

    math::Point
    center(
            AtomOPs const &atoms) {
        auto center = math::Point();
        for (auto const &a : atoms) {
            center += a->get_coords();
        }
        return center / float(atoms.size());
    }

    void
    Residue::_build_beads() {
        if (res_type_->get_set_type() == SetType::RNA) { _build_beads_RNA(); }
        else if (res_type_->get_set_type() == SetType::PROTEIN) {
            beads_.push_back(util::Bead(get_coords("CA"), util::BeadType::CALPHA));
        } else {
            beads_.push_back(util::Bead(get_center(), util::BeadType::MCENTER));
        }
    }

    void
    Residue::_build_beads_RNA() {
        std::vector<Atom const *> phos_atoms, sugar_atoms, base_atoms;
        int i = -1;
        for (auto const &a : atoms_) {
            i++;
            if (i < 3) { phos_atoms.push_back(&a); }
            else if (i < 12) { sugar_atoms.push_back(&a); }
            else { base_atoms.push_back(&a); }
        }

        auto get_center = [&](
                std::vector<Atom const *> const &atom_ptrs) {
            auto center = math::Point();
            for (auto a : atom_ptrs) {
                center = center + a->get_coords();
            }
            center = center / float(atom_ptrs.size());
            //std::cout << atom_ptrs.size() << std::endl;

            return center;
        };

        if (phos_atoms.size() > 0) {
            beads_.push_back(util::Bead(get_center(phos_atoms), util::BeadType::PHOS));
        }
        if (sugar_atoms.size() > 0) {
            beads_.push_back(util::Bead(get_center(sugar_atoms), util::BeadType::SUGAR));
        }
        if (base_atoms.size() > 0) {
            beads_.push_back(util::Bead(get_center(base_atoms), util::BeadType::BASE));
        }
    }

    String
    Residue::get_str() const {
        std::stringstream ss;
        ss << res_type_->get_name() << "," << name_ << "," << num_ << "," << chain_id_ << "," << i_code_ << ",";
        for (auto const &a : atoms_) {
            ss << a.get_str() + ",";
        }
        return ss.str();
    }

    json::JSON
    Residue::get_json() const {
        auto atom_array = json::Array();
        for (auto const &a : atoms_) { atom_array.append(a.get_json()); }
        return json::JSON{
                "res_type", res_type_->get_name(),
                "name", name_,
                "num", num_,
                "chain_id", chain_id_,
                "i_code", i_code_,
                "atoms", atom_array};
    }


    String
    Residue::get_pdb_str(
            int &acount,
            int rnum,
            char chain_id) const {

        auto s = String();
        for (auto const &a : atoms_) {
            char buffer[200];
            std::sprintf(buffer, "%-6s%5d %-4s%1s%-4c%1c%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s\n", "ATOM",
                         acount, a.get_str_name().c_str(), "", name_, chain_id, rnum, "", a.get_x(), a.get_y(),
                         a.get_z(), 1.00, 0.00, "", "");
            s += String(buffer);
            acount++;
        }
        return s;
    }

    String
    Residue::get_bead_pdb_str(
            int &acount,
            int rnum,
            char chain_id) const {
        auto s = String();
        auto bead_name = String("C");
        for (auto const &b : beads_) {
            char buffer[200];
            auto c = b.get_center();
            if (b.get_type() == util::BeadType::BASE) { bead_name = "C"; }
            else if (b.get_type() == util::BeadType::SUGAR) { bead_name = "N"; }
            else if (b.get_type() == util::BeadType::PHOS) { bead_name = "P"; }
            std::sprintf(buffer, "%-6s%5d %-4s%1s%-4c%1c%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s\n", "ATOM",
                         acount, bead_name.c_str(), "", name_, chain_id, rnum, "", c.get_x(), c.get_y(), c.get_z(),
                         1.00, 0.00, "", "");
            s += String(buffer);
            acount++;

        }
        return s;
    }


    void
    Residue::write_pdb(
            String const fname) {
        std::ofstream out;
        out.open(fname.c_str());
        out << get_pdb_str(1) << std::endl;
        out.close();
    }

}
