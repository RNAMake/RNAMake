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
#include "structure/residue.h"

Point
calc_center(AtomOPs const & atoms) {
    assert(atoms.size() > 0);
    Point center(0, 0, 0);
    auto total = 0;
    for (auto const & a : atoms) {
        if (a == nullptr) {
            continue;
        }
        center += a->coords();
        total += 1;
    }
    return center / total;
}

Residue::Residue(
        AtomOPs const & atoms,
        ResidueTypeOP const & rtype,
        String const & name,
        int const & num,
        String const & chain_id,
        String const & i_code) :
        primitives::Residue(name, num, chain_id, i_code),
        rtype_(rtype),
        beads_(Beads()) {

    setup_atoms(atoms);
}

Residue::Residue(
        Residue const & r,
        int new_uuid,
        int build_beads):
        primitives::Residue(r.name_, r.num_, r.chain_id_, r.i_code_, r.uuid_) {

    rtype_ = r.rtype_;
    atoms_ = AtomOPs(r.atoms_.size());
    int i = -1;
    for (auto const & a : r.atoms_) {
        i++;
        if (a == nullptr) { continue; }
        atoms_[i] = std::make_shared<Atom>(*a);
    }

    if(new_uuid) { uuid_ = Uuid(); }
    if(build_beads && r.beads_.size() > 0) {
        beads_ = _get_beads();
    }
}

Residue::Residue(
        Residue const & r,
        Uuid const & given_uuid,
        int build_beads) :
        primitives::Residue(r.name_, r.num_, r.chain_id_, r.i_code_, given_uuid) {

    rtype_ = r.rtype_;
    atoms_ = AtomOPs(r.atoms_.size());
    int i = -1;
    for (auto const & a : r.atoms_) {
        i++;
        if (a == nullptr) { continue; }
        atoms_[i] = std::make_shared<Atom>(*a);
    }

    if (build_beads && r.beads_.size() > 0) {
        beads_ = _get_beads();
    }
}

Residue::Residue(
        String const & s,
        ResidueTypeSet const & rts):
        primitives::Residue() {

    Strings spl = split_str_by_delimiter(s, ",");
    rtype_    = rts.get_type(spl[0]);
    name_     = spl[1];
    num_      = std::stoi(spl[2]);
    chain_id_ = spl[3];
    i_code_   = spl[4];
    atoms_    = AtomOPs();
    auto atoms = AtomOPs();
    int i = 5;
    while (i < spl.size()) {
        if (spl[i].length() == 1) {
            atoms.push_back(nullptr);
        } else {
            atoms.push_back(std::make_shared<Atom>(spl[i]));
        }
        i++;
    }
    setup_atoms(atoms);
}


//setup helpers ////////////////////////////////////////////////////////////////////////////////////

void
Residue::setup_atoms(
        AtomOPs const & atoms) {
    atoms_ = AtomOPs(rtype_->size());
    int count = 0;
    for (auto const & a : atoms) {
        if (a == nullptr) { continue; }
        // does this atom belong in this residue
        if(!rtype_->is_valid_atom(a->name())) { continue; }

        auto name_change = rtype_->get_correct_atom_name(*a);
        //check for misnamed atoms
        if(name_change.length() != 0) {
            auto new_a = std::make_shared<Atom>(name_change, a->coords());
            int pos = rtype_->atom_index(new_a->name());
            if (pos == -1) { continue; }
            atoms_[pos] = new_a;
        }
        else {
            int pos = rtype_->atom_index(a->name());
            if (pos == -1) { continue; }
            atoms_[pos] = a;
        }
        count++;
    }

    if (count == 0) { throw ResidueException("Residue has zero atoms something wrong"); }
}

Beads
Residue::_get_beads() const {
    AtomOPs phos_atoms, sugar_atoms, base_atoms;
    int i = -1;
    for (auto const & a : atoms_) {
        i++;
        if (a == nullptr) { continue; }
        if (i < 3) { phos_atoms.push_back(a); }
        else if (i < 12) { sugar_atoms.push_back(a); }
        else { base_atoms.push_back(a); }
    }
    auto beads = Beads();
    if (phos_atoms.size() > 0)  { beads.push_back(Bead(calc_center(phos_atoms), BeadType::PHOS)); }
    if (sugar_atoms.size() > 0) { beads.push_back(Bead(calc_center(sugar_atoms), BeadType::SUGAR)); }
    if (base_atoms.size() > 0)  { beads.push_back(Bead(calc_center(base_atoms), BeadType::BASE)); }
    return beads;
}

//public methods ///////////////////////////////////////////////////////////////////////////////////

String
Residue::to_str() const {
    std::stringstream ss;
    ss << rtype_->name() << "," << name_ << "," << num_ << "," << chain_id_ << "," << i_code_ << ",";
    for (auto const & a : atoms_) {
        if (a == nullptr) { ss << "N,"; }
        else { ss << a->to_str() + ","; }
    }
    return ss.str();
}

String
Residue::to_pdb_str(
        int & acount,
        int rnum,
        String const & chain_id) const {

    int num = num_;
    if (rnum != -1) {
        num = rnum;
    }
    String cid = chain_id_;
    if (chain_id != "") {
        cid = chain_id;
    }

    String s;
    for (auto const & a : atoms_) {
        if (a == nullptr) { continue; }
        char buffer[200];
        std::sprintf(buffer, "%-6s%5d %-4s%1s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s\n", "ATOM", acount,
                     a->name().c_str(), "", short_name().c_str(), cid.c_str(), num, "", a->coords()[0], a->coords()[1],
                     a->coords()[2], 1.00, 0.00, "", "");
        s += String(buffer);
        acount++;
    }
    return s;
}

void
Residue::to_pdb(
        String const fname) {
    std::ofstream out;
    out.open(fname.c_str());
    int i = 1;
    String s = to_pdb_str(i);
    out << s << std::endl;
    out.close();
}

state::ResidueOP
Residue::get_state() {
    return std::make_shared<state::Residue>(name_, num_, chain_id_, i_code_, beads_, uuid_);
}






















