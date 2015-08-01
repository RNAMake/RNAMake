//
//  chain.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <memory.h>

#include "base/string.h"
#include "secondary_structure/chain.h"

namespace sstruct {

Chain
str_to_chain(String const & s) {
    Strings spl = split_str_by_delimiter(s, ";");
    ResidueOPs res;
    for(auto const & r_str : spl) {
        if(r_str.length() < 3) { continue; }
        auto r = std::make_shared<Residue>(str_to_residue(r_str));
        res.push_back(r);
    }
    return Chain(res);
}

}