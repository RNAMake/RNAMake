//
//  residue_type.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

//RNAMake Headers
#include "structure/residue_type.h"

ResidueType::ResidueType(
        String const & name,
        StringIntMap const & atom_map,
        SetType const & set_type):
        name_(name),
        atom_map_(atom_map),
        set_type_(set_type) {

    alt_names_ = Strings();
    if (set_type_ == SetType::RNA) {
        alt_names_.push_back(name_.substr(0, 1));
        alt_names_.push_back("r" + name_.substr(0, 1));
        alt_names_.push_back("D" + name_.substr(0, 1));
        extend_res_specific_altnames();
    }

    atom_alt_names_ = StringStringMap();
    atom_alt_names_["O1P"] = "OP1";
    atom_alt_names_["O2P"] = "OP2";
}

bool
ResidueType::is_valid_atom(String const & name) const{
    if(atom_map_.find(name) == atom_map_.end()) { return false; }
    else { return true; }
}

int
ResidueType::atom_index(String const & name) const {
    if(atom_map_.find(name) != atom_map_.end()) { return atom_map_.at(name); }
    throw ResidueTypeException("atom: " + name + " does not exist in ResidueType");
}

String
ResidueType::get_correct_atom_name(Atom const & a) const  {

    if(atom_alt_names_.find(a.name()) == atom_alt_names_.end()) { return String(""); }
    return atom_alt_names_.at(a.name());
}

int
ResidueType::match_name(
        String const & name) const {
    if (name_ == name ) { return 1; }
    for (auto const & s : alt_names_) {
        if (s == name) { return 1; }
    }
    return 0;
}

void
ResidueType::extend_res_specific_altnames() {
    // There has to be a better way to do this
    Strings alt_names;
    if (name_ == "GUA") {
        alt_names = split_str_by_delimiter("MIA GDP GTP M2G 1MG 7MG G7M QUO I YG", " ");
    } else if (name_ == "ADE") {
        alt_names = split_str_by_delimiter("A23 3DA 1MA 12A AET 2MA", " ");
    } else if (name_ == "URA") {
        alt_names = split_str_by_delimiter("PSU H2U 5MU 4SU 5BU 5MC U3H 2MU 70U BRU DT", " ");
    } else if (name_ == "CYT") {
        alt_names = split_str_by_delimiter("CBR CCC", " ");
    } else {

    }

    for (auto const & name : alt_names) {
        alt_names_.push_back(name);
    }


}






