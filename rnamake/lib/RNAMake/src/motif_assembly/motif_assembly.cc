//
//  motif_assembly.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_assembly.h"

void
MotifAssembly::setup_options() {
    
    options_ = Options();
    options_.add_option(Option("error_on_add", 1));
    
}

MotifTreeNodeOP
MotifAssembly::add_motif(
    MotifOP const & m,
    MotifTreeNodeOP parent,
    String end_name,
    String parent_name) {
    
    
    
    MotifTreeNodeOP mt_node = mt_.add_motif(m);
    
    if(mt_node == nullptr && option<int>("error_on_add") == 1) {
        throw "Cannot add motif: " + m->name() + " to assembly";
    }
    
    return mt_node;
    
}