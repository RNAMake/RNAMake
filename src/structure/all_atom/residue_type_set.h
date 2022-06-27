//
//  residue_type_set.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__residue_type_set__
#define __RNAMake__residue_type_set__

#include <cstdio>
#include <structure/all_atom/residue_type.h>

namespace structure::all_atom {

class ResidueTypeSet {
public:
  ResidueTypeSet();

  ~ResidueTypeSet() {}

  const ResidueType &get_residue_type(String const &) const;

  bool contains_residue_type(String const &) const;

private:
  String _get_rtype_name(String const &);

  StringIntMap _get_atom_map_from_file(String const &);

  void _read_rtypes_from_dir(String const &, SetType const &);

  Strings _get_extra_resnames_for_specific_res(String const &);

private:
  std::vector<ResidueType> _residue_types = {};
};

} // namespace structure::all_atom

#endif /* defined(__RNAMake__residue_type_set__) */