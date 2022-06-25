//
//  residue_type.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

// RNAMake Headers
#include <structure/residue_type.h>

ResidueType::ResidueType(String const &name, StringIntMap const &atom_name_map,
                         SetType set_type, Strings const &extra_alt_names)
    : _name(name), _atom_name_map(atom_name_map), _alt_names(extra_alt_names),
      _set_type(set_type) {

  if (_set_type == SetType::RNA) {
    _alt_names.push_back(_name.substr(0, 1));
    _alt_names.push_back("r" + _name.substr(0, 1));
    _alt_names.push_back("D" + _name.substr(0, 1));
  }
}

bool ResidueType::is_valid_atom_name(String const &atom_name) const {
  if (_atom_name_map.find(atom_name) != _atom_name_map.end()) {
    return true;
  } else {
    return false;
  }
}

Index ResidueType::get_atom_index(String const &name) const {
  expects<ResidueTypeException>(is_valid_atom_name(name),
                                "must supply valid atom name for residue type");
  return _atom_name_map.find(name)->second;
}

String ResidueType::get_atom_name_at_pos(Index index) const {
  for (auto const &kv : _atom_name_map) {
    if (kv.second == index) {
      return kv.first;
    }
  }
  throw ResidueTypeException(
      "residue type: " + get_name() +
      " does not have an index: " + std::to_string(index));
}

bool ResidueType::is_valid_residue_name(String const &res_name) const {
  if (res_name == _name) {
    return true;
  }
  for (auto const &alt_name : _alt_names) {
    if (res_name == alt_name) {
      return true;
    }
  }
  return false;
}

ResidueTypeCOP get_new_residue_type(String const &res_name,
                                    Strings const &atom_names) {
  expects<ResidueTypeException>(
      atom_names.size() != 0,
      "must have more than one atom name to create new atom type");

  auto atom_name_map = StringIntMap();
  int i = 0;
  for (auto const &name : atom_names) {
    atom_name_map[name] = i;
    i++;
  }
  auto extra_alt_names = Strings();

  return std::make_shared<ResidueType>(res_name, atom_name_map,
                                       SetType::UNKNOWN, extra_alt_names);
}
