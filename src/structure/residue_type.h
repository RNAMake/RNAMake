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
#include "base/types.h"
#include "structure/atom.h"

namespace structure {

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
    ResidueType() {}

    ResidueType(
            String const &,
            StringIntMap const &,
            SetType const & set_type);

    ~ResidueType() {}

public:
    String
    get_correct_atom_name(
            Atom const &) const;

    int
    match_name(
            String const &) const;

    inline
    String
    short_name() const { return name_.substr(0, 1); }

    inline
    int
    atom_pos_by_name(
            String const & aname) const {

        StringIntMap::const_iterator iter(atom_map_.find(aname));
        if (iter != atom_map_.end()) {
            return iter->second;
        } else {
            return -1;
        }
    }

    inline
    int
    size() { return (int) atom_map_.size(); }

public: //getters

    inline
    String
    const &
    name() const { return name_; }

    inline
    SetType
    set_type() const { return set_type_; }


private:

    void
    extend_res_specific_altnames();

private:
    String name_;
    StringIntMap atom_map_;
    Strings alt_names_;
    StringStringMap atom_alt_names_;
    SetType set_type_;

};

typedef std::vector<ResidueType> ResidueTypes;

}

#endif /* defined(__RNAMake__residue_type__) */
