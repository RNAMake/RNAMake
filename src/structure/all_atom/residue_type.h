//
//  residue_type.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__residue_type__
#define __RNAMake__residue_type__

#include <cstdio>

// RNAMake Headers
//#include <base/assertions.h>
#include <base/types.hpp>

namespace structure::all_atom {

/// @brief - exception for residue type
class ResidueTypeException : public std::runtime_error {
public:
  /**
   * Standard constructor for ResidueException
   * @param   message   Error message for residue type
   */
  ResidueTypeException(String const &message) : std::runtime_error(message) {}
};

/// @brief - enum of residue types
enum class SetType { RNA, PROTEIN, UNKNOWN };

class ResidueType {
  /// @brief - constructor
public:
  ResidueType(String const &, StringIntMap const &, SetType, Strings const &);

  ~ResidueType() = default;

  /// @brief - checks if the atom name is valid/something real
  bool is_valid_atom_name(String const &) const;

  /// @brief - gets the index of an atom from a given string of a residue
  Index get_atom_index(String const &) const;

  /// @brief - gets the name of an atom at a specified position/index
  String get_atom_name_at_pos(Index) const;

  /// @brief - checks if the residue name is valid
  bool is_valid_residue_name(String const &) const;

public: // getters
  /// @brief - gets the name of the residue type
  inline String const &get_name() const { return _name; }

  /// @brief - gets the shortname of the residue type
  inline char get_short_name() const { return _name[0]; }

  /// @brief - gets set type
  inline SetType get_set_type() const { return _set_type; }

  /// @brief - gets the number of atoms in the residue
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

} // namespace structure::all_atom
#endif /* defined(__RNAMake__residue_type__) */