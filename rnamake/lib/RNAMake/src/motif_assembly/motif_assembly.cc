//
//  motif_assembly.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_assembly.h"
#include "resources/library_manager.h"

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
    String parent_name,
    int end_index) {
    
    int parent_end_index = -1;
    
    if(end_name.length() > 0) { end_index = m->end_index(m->get_basepair_by_name(end_name)); }
    if(parent_name.length() > 0) {
        if(parent == nullptr) {
            parent = mt_.last_node();
        }
        
        parent_end_index = parent->motif()->end_index(parent->motif()->get_basepair_by_name(end_name));
    }

    MotifTreeNodeOP mt_node = mt_.add_motif(m, parent, end_index, 0, parent_end_index);

    if(mt_node == nullptr && option<int>("error_on_add") == 1) {
        throw "Cannot add motif: " + m->name() + " to assembly";
    }
    
    return mt_node;
    
}

MotifTreeNodeOP
MotifAssembly::add_motif(
    String const & name,
    MotifTreeNodeOP parent,
    String end_name,
    String parent_name,
    int end_index) {
    
    MotifOP m = LibraryManager::getInstance().get_motif(name);
    return add_motif(m, parent, end_name, parent_name, end_index );
    
}