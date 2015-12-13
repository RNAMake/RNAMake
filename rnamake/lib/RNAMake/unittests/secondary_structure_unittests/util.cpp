//
//  util.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/10/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "util.h"
#include "util/random_number_generator.h"
#include "secondary_structure/util.h"


namespace unittests {
namespace sstruct_unittests  {

void
fill_basepairs_in_ss(sstruct::PoseOP & ss) {
    auto rng = RandomNumberGenerator();
    auto pairs = Strings{"AU", "UA", "GC", "CG"};
    
    for(auto & bp : ss->basepairs()) {
        if(bp->res1()->name() != "N" || bp->res2()->name() != "N") { continue; }
        auto pos = rng.randrange(4);
        auto p = pairs[pos];
        String name1, name2;
        name1.push_back(p[0]); name2.push_back(p[1]);
        bp->res1()->name(name1);
        bp->res2()->name(name2);
    }
    
    for(auto & m : ss->motifs()) {
        auto end_ids = Strings();
        for(auto const & end : m->ends()) {
            end_ids.push_back(sstruct::assign_end_id(m, end));
        }
        m->end_ids(end_ids);
    }
}

}
}
