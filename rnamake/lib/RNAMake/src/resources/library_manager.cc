//
//  library_manager.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "util/settings.h"
#include "motif/motif_type.h"
#include "motif/motif_tree.h"
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
    throw "cannot find motif " + name + "\n";
    
}

void
LibraryManager::add_motif(String const & path) {
    String name = filename(path);
    //extra_motifs_[name] =
    
}