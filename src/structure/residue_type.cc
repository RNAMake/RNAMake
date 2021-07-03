//
//  residue_type.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

//RNAMake Headers
#include <structure/residue_type.h>


ResidueType::ResidueType(
        String const & name,
        StringIntMap const & atom_name_map,
        SetType set_type,
        Strings const & extra_alt_names):
        name_(name),
        atom_name_map_(atom_name_map),
        alt_names_(extra_alt_names),
        set_type_(set_type) {

    if (set_type_ == SetType::RNA) {
        alt_names_.push_back(name_.substr(0, 1));
        alt_names_.push_back("r" + name_.substr(0, 1));
        alt_names_.push_back("D" + name_.substr(0, 1));
    }
}

bool
ResidueType::is_valid_atom_name(
        String const & atom_name)  const {
    if(atom_name_map_.find(atom_name) != atom_name_map_.end()) { return true; }
    else                                                       { return false; }
}

Index
ResidueType::get_atom_index(
        String const & name) const {
    expects<ResidueTypeException>(
            is_valid_atom_name(name),
            "must supply valid atom name for residue type");
    return atom_name_map_.find(name)->second;
}

String
ResidueType::get_atom_name_at_pos(
        Index index) const {
    for(auto const & kv : atom_name_map_) {
        if(kv.second == index) { return kv.first; }
    }
    throw ResidueTypeException("residue type: " + get_name() +" does not have an index: " + std::to_string(index));
}

bool
ResidueType::is_valid_residue_name(
        String const & res_name) const {
    if(res_name == name_) { return true; }
    for(auto const & alt_name : alt_names_) {
        if(res_name == alt_name) { return true; }
    }
    return false;
}

ResidueTypeCOP
get_new_residue_type(
        String const & res_name,
        Strings const & atom_names) {
    expects<ResidueTypeException>(
            atom_names.size() != 0,
            "must have more than one atom name to create new atom type");

    auto atom_name_map = StringIntMap();
    int i = 0;
    for(auto const & name : atom_names) {
        atom_name_map[name] = i;
        i++;
    }
    auto extra_alt_names = Strings();

    return std::make_shared<ResidueType>(res_name, atom_name_map, SetType::UNKNOWN, extra_alt_names);
}




