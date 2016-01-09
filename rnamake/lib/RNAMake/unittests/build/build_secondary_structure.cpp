//
//  build_secondary_structure.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "util/random_number_generator.h"
#include "secondary_structure/secondary_structure_factory.h"
#include "build_secondary_structure.h"

namespace unittests {
    
sstruct::PoseOP
BuildSecondaryStructure::build_helix(
    int size,
    int opt) {
    auto rng = RandomNumberGenerator();
    auto pairs = Strings{"AU", "UA", "GC", "CG"};
    String s1, s2, ss1, ss2;
    
    for(int i = 0; i < size; i++) {
        auto pos = rng.randrange(4);
        auto p = pairs[pos];
        
        s1 += p[0];
        s2 = p[1] + s2 ;
        ss1 += "(";
        ss2 += ")";
        
    }
    
    auto ssf = sstruct::SecondaryStructureFactory();
    return ssf.pose(s1+"&"+s2, ss1+"&"+ss2);
}

sstruct::PoseOP
BuildSecondaryStructure::build_hairpin(int size, int opt) {
    auto rng = RandomNumberGenerator();
    auto pairs = Strings{"AU", "UA", "GC", "CG"};
    if(opt) {
        pairs = Strings{"NN"};
    }
    
    String s1, s2, ss1, ss2;
    
    for(int i = 0; i < size; i++) {
        auto pos = rng.randrange(pairs.size());
        auto p = pairs[pos];
        
        s1 += p[0];
        s2 = p[1] + s2 ;
        ss1 += "(";
        ss2 += ")";
        
    }
    
    auto ssf = sstruct::SecondaryStructureFactory();
    return ssf.pose(s1+"AAAA"+s2, ss1+"...."+ss2);
}
    
    
}