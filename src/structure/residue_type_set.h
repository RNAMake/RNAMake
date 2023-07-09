//
//  residue_type_set.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__residue_type_set__
#define __RNAMake__residue_type_set__

#include <stdio.h>
#include "residue_type.h"

namespace structure {

class ResidueTypeSet {
public:
    ResidueTypeSet();

    ~ResidueTypeSet() {}

    ResidueType const &
    get_rtype_by_resname(
            String const &) const;

    bool
    contains_rtype(
            String const &) const;

private:
    String
    _get_rtype_name(
            String const &);

    StringIntMap
    _get_atom_map_from_file(
            String const &);

    void
    _read_rtypes_from_dir(
            String const &,
            SetType const &);


private:
    ResidueTypes residue_types_;

};

}

#endif /* defined(__RNAMake__residue_type_set__) */
