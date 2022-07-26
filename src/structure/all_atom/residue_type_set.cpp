//
//  residue_type_set.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <algorithm>
#include <dirent.h>
#include <fstream>

// RNAMake Headers
#include <base/paths.hpp>
#include <base/string.hpp>
#include <structure/all_atom/residue_type_set.h>

namespace structure::all_atom {

ResidueTypeSet::ResidueTypeSet() {
  auto base_path = base::path::resources_path() + "/residue_types/";
  _read_rtypes_from_dir(base_path + "RNA/", SetType::RNA);
  _read_rtypes_from_dir(base_path + "PROTEIN/", SetType::PROTEIN);
}

void ResidueTypeSet::_read_rtypes_from_dir(String const &path,
                                           SetType const &set_type) {
  DIR *pDIR;
  struct dirent *entry;
  pDIR = opendir(path.c_str());
  while ((entry = readdir(pDIR)) != nullptr) {
    auto fname = String(entry->d_name);
    if (fname.length() < 4) {
      continue;
    }
    auto name = _get_rtype_name(path + fname);
    auto atom_map = _get_atom_map_from_file(path + fname);
    auto alt_names = _get_extra_resnames_for_specific_res(name);
    _residue_types.emplace_back(name, atom_map, set_type, alt_names);
  }
  closedir(pDIR);
  delete entry;
}

String ResidueTypeSet::_get_rtype_name(String const &fname) {
  auto name_spl = base::string::split(fname, "/");
  auto type_file_name = name_spl.back();
  auto spl_2 = base::string::split(type_file_name, ".");
  return spl_2[0];
}

StringIntMap ResidueTypeSet::_get_atom_map_from_file(String const &fname) {

  auto line = String();
  std::ifstream input;
  input.open(fname);
  getline(input, line);
  input.close();
  auto atom_names = base::string::split(line, " ");
  auto i = 0;
  auto atom_map = StringIntMap();
  for (auto const &name : atom_names) {
    atom_map[name] = i;
    i++;
  }
  return atom_map;
}

const ResidueType &
ResidueTypeSet::get_residue_type(String const &resname) const {
  std::string res = "blank";
  for (auto const &restype : _residue_types) {
    if (resname[0] == ':') {
      // TODO remove these string manipulation, string should be passed as
      // expected.
      // res = replace_char(resname, ':', ' ');
      // res.erase(0, 1);
      if (restype.is_valid_residue_name(res)) {
        return restype;
      }
    } else {
      res = resname;
      if (restype.is_valid_residue_name(resname)) {
        return restype;
      }
    }
  }

  throw ResidueTypeException("cannot find residue with name: " + res);
}

bool ResidueTypeSet::contains_residue_type(String const &resname) const {
  for (auto const &restype : _residue_types) {
    if (restype.is_valid_residue_name(resname)) {
      return true;
    }
  }
  return false;
}

Strings
ResidueTypeSet::_get_extra_resnames_for_specific_res(String const &res_name) {
  auto alt_names = Strings();
  if (res_name == "GUA") {
    alt_names =
        base::string::split("MIA GDP GTP M2G 1MG 7MG G7M QUO I YG", " ");
  } else if (res_name == "ADE") {
    alt_names = base::string::split("A23 3DA 1MA 12A AET 2MA", " ");
  } else if (res_name == "URA") {
    alt_names =
        base::string::split("PSU H2U 5MU 4SU 5BU 5MC U3H 2MU 70U BRU DT", " ");
  } else if (res_name == "CYT") {
    alt_names = base::string::split("CBR CCC", " ");
  }
  return alt_names;
}

} // namespace structure::all_atom