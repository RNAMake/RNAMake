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

//RNAMake Headers
#include <base/types.h>
#include <base/assertions.h>
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
    ResidueTypeException(String const & message) :
            std::runtime_error(message) {}
};

enum class SetType {
    RNA, PROTEIN, UNKNOWN
};

class ResidueType {
public:
    ResidueType(
            String const &,
            StringIntMap const &,
            SetType,
            Strings const &);

    ~ResidueType() {}

public:
    bool
    is_valid_atom_name(
            String const &) const;

    Index
    get_atom_index(
            String const &) const;

    String
    get_atom_name_at_pos(
            Index) const;

    bool
    is_valid_residue_name(
            String const &) const;

public: //getters

    inline
    String
    const &
    get_name() const { return name_; }

    inline
    char
    get_short_name() const { return name_[0]; }

    inline
    SetType
    get_set_type() const { return set_type_; }

    inline
    size_t
    get_num_atoms() const { return atom_name_map_.size(); }
private:


private:
    String name_;
    StringIntMap atom_name_map_;
    Strings alt_names_;
    SetType set_type_;

};

typedef std::shared_ptr<ResidueType>       ResidueTypeOP;
typedef std::shared_ptr<ResidueType const> ResidueTypeCOP;
typedef std::vector<ResidueTypeOP>         ResidueTypeOPs;

ResidueTypeCOP
get_new_residue_type(
        String const & res_name,
        Strings const & atom_names);


#endif /* defined(__RNAMake__residue_type__) */