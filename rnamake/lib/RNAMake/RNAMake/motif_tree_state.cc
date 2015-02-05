//
//  motif_tree_state.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_state.h"
#include "types.h"
#include "FileIO.h"

NameElements
parse_db_name(
    String const & s) {
    Strings spl = split_str_by_delimiter(s, "-");
    return NameElements(spl[0], std::stoi(spl[1]), std::stoi(spl[2]),
                        std::stoi(spl[3]), std::stoi(spl[4]), std::stoi(spl[5]), std::stoi(spl[6]));
}