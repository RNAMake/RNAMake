//
//  library_manager.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif/motif_type.h"
#include "resources/motif_library.h"
#include "resources/library_manager.h"


LibraryManager::LibraryManager() {
    mlibs_ = std::map<String, MotifLibraryOP>();
    mlibs_["HELIX"   ] = MotifLibraryOP(new MotifLibrary(HELIX));
    mlibs_["NWAY"    ] = MotifLibraryOP(new MotifLibrary(NWAY));
    mlibs_["TWOWAY"  ] = MotifLibraryOP(new MotifLibrary(TWOWAY));
    mlibs_["TCONTACT"] = MotifLibraryOP(new MotifLibrary(TCONTACT));
}


MotifOP
LibraryManager::get_motif(String const & name) {
    for(auto const & kv : mlibs_) {
        if(kv.second->contains_motif(name)) {
            return kv.second->get_motif(name);
        }
    }
    throw "cannot find motif " + name + "\n";
    
}