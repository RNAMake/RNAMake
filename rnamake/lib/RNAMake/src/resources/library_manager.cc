//
//  library_manager.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//


//RNAMake Headers
#include "util/settings.h"
#include "motif/motif_type.h"
#include "motif/motif_tree.h"
#include "motif/motif.h"
#include "resources/motif_library.h"
#include "resources/library_manager.h"


LibraryManager::LibraryManager() {
    mlibs_ = std::map<String, MotifLibraryOP>();
    extra_motifs_ = std::map<String, MotifOP>();
    mlibs_["HELIX"   ] = MotifLibraryOP(new MotifLibrary(HELIX));
    mlibs_["NWAY"    ] = MotifLibraryOP(new MotifLibrary(NWAY));
    mlibs_["TWOWAY"  ] = MotifLibraryOP(new MotifLibrary(TWOWAY));
    mlibs_["TCONTACT"] = MotifLibraryOP(new MotifLibrary(TCONTACT));
    
    String path = resources_path() + "motif_libraries/bp_steps.db";
    MotifLibraryOP mlib (new MotifLibrary(path, HELIX));
    mlibs_["BP_STEPS"] = mlib;
    
    mt_  = MotifTree();
    mt2_ = MotifTree();
    mt_.add_motif(get_motif("HELIX.IDEAL.6"), nullptr, 1, -1, 0);
    mt_.increase_level();
    
}


MotifOP
LibraryManager::get_motif(
    String const & name,
    int const & end_index,
    String const & end_name) {
    for(auto const & kv : mlibs_) {
        if(kv.second->contains_motif(name)) {
            return kv.second->get_motif(name);
        }
    }
    
    if(extra_motifs_.find(name) != extra_motifs_.end()) {
        if(end_index == -1 && end_name.length() == 0) {
            return extra_motifs_[name];
        }
    }
    
    throw "cannot find motif " + name + "\n";
    
}

void
LibraryManager::add_motif(String const & path) {
    String name = filename(path);
    extra_motifs_[name] = MotifOP(new Motif(path));
    
}

MotifOP
LibraryManager::_prep_extra_motif_for_asssembly(
    MotifOP const & m,
    int end_index,
    String const & end_name) {
    
    if(end_index == -1 && end_name.length() == 0) {
        throw "cannot call _prep_extra_motif_for_asssembly without end_name or end_index";
    }
    
    if(end_name.length() != 0) {
        BasepairOP end = m->get_basepair_by_name(end_name);
        end_index = m->end_index(end);
    }
    
    MotifTreeNodeOP node = mt_.add_motif(m, nullptr, end_index);
    if(node == nullptr) {
        throw "cannot prep motif for assembly, motif: " + m->name();
    }
    
    MotifOP h_m = get_motif("HELIX.IDEAL.3");
    for(auto const & e : node->available_ends()) {
        MotifTreeNodeOP n = mt_.add_motif(h_m, node, 1, -1, 0);
        if(n == nullptr) {
            e->flip();
        }
        else {
           // mt_.remove_node(mt_.last_node());
        }
    }
     
    return m;
    
    
}