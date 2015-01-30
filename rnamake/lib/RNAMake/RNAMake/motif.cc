//
//  motif.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif.h"

Motif::Motif(
    String const & s,
    ResidueTypeSet const & rts):
    beads_(Beads()),
    score_(0),
    basepairs_(Basepairs()),
    ends_(Basepairs()),
    mdir_(String()),
    name_(String()),
    cached_rotations_(Matrices())
{
    
    Strings spl = split_str_by_delimiter(s, "&");
    mdir_ = spl[0];
    name_ = spl[1];
    score_ = std::stof(spl[2]);
    mtype_ = static_cast<MotifType>(std::stoi(spl[3]));
    structure_ = str_to_structure(spl[4], rts);
    Strings basepair_str = split_str_by_delimiter(spl[5], "@");
    for (auto const & bp_str : basepair_str) {
        Strings bp_spl = split_str_by_delimiter(bp_str, ",");
        Strings res_spl = split_str_by_delimiter(bp_spl[0], "-");
        String res1_id = res_spl[0].substr(0,1);
        int res1_num = std::stoi(res_spl[0].substr(1));
        Residue res1 = structure_.get_residue(res1_num, res1_id, "");
    }
    
}
