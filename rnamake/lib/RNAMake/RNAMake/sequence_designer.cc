//
//  sequence_designer.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "sequence_designer.h"

String
SequenceDesigner::design(
    String const & sequence,
    String const & structure) {
    
    
    ss_tree_ = SecondaryStructureTree(structure, sequence);
    SecondaryStructureNodeOPs designable_bps = ss_tree_.get_designable_bps();
        
    return "";
    
}


void
mutate_basepair(
    SecondaryStructureNodeOP & bp_node) {
    
    
}






