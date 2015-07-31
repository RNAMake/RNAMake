//
//  residue.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "residue.h"

namespace sstruct {

Residue
str_to_residue(String const & s) {
    Strings spl = split_str_by_delimiter(s, ",");
    return Residue(spl[0], spl[1], std::stoi(spl[2]), spl[3], Uuid(), spl[4]);
}

}