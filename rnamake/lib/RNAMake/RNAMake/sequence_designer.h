//
//  sequence_designer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sequence_designer__
#define __RNAMake__sequence_designer__

#include <stdio.h>
#include "random_number_generator.h"
#include "secondary_structure_tree.h"
#include "vienna.h"
#include "FileIO.h"

class SequenceDesigner {
public:
    SequenceDesigner():
    v_ (Vienna() ) {
        
        rng_  = RandomNumberGenerator();
        basepairs_ = split_str_by_delimiter("GC,CG,AU,UA", ",");
    }
    
    ~SequenceDesigner() {}
    
public:
    
    String
    design(
        String const &,
        String const &);
    
    void
    mutate_basepair(
        SecondaryStructureNodeOP & bp_node);
    
private:
    RandomNumberGenerator rng_;
    SecondaryStructureTree ss_tree_;
    Vienna v_;
    Strings basepairs_;
    String last_bp_type_;
};

#endif /* defined(__RNAMake__sequence_designer__) */
