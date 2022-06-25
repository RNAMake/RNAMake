//
//  residue_type.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__residue_type__
#define __RNAMake__residue_type__

#include <stdio.h>
#include <vector>

// RNAMake Headers
//#include <base/assertions.h>
#include <base/types.hpp>
#include <structure/atom.h>

/*
 * Exception for residue type
 */
class ResidueTypeException : public std::runtime_error {
public:
  /**
   * Standard constructor for ResidueException
   * @param   message   Error message for residue type
   */
  ResidueTypeException(String const &message) : std::runtime_error(message) {}
};

enum class SetType { RNA, PROTEIN, UNKNOWN };

class ResidueType {
public:
  ResidueType(String const &, StringIntMap const &, SetType, Strings const &);

  ~ResidueType() {}

public:
  bool is_valid_atom_name(String const &) const;

  Index get_atom_index(String const &) const;

  String get_atom_name_at_pos(Index) const;

  bool is_valid_residue_name(String const &) const;

public: // getters
  inline String const &get_name() const { return _name; }

  inline char get_short_name() const { return _name[0]; }

  inline SetType get_set_type() const { return _set_type; }

  inline size_t get_num_atoms() const { return _atom_name_map.size(); }

private:
private:
  String _name;
  StringIntMap _atom_name_map;
  Strings _alt_names;
  SetType _set_type;
};

typedef std::shared_ptr<ResidueType> ResidueTypeOP;
typedef std::shared_ptr<ResidueType const> ResidueTypeCOP;
typedef std::vector<ResidueTypeOP> ResidueTypeOPs;

ResidueTypeCOP get_new_residue_type(String const &res_name,
                                    Strings const &atom_names);

#endif /* defined(__RNAMake__residue_type__) */