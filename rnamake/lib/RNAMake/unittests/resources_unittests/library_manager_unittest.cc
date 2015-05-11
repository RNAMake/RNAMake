//
//  library_manager_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "library_manager_unittest.h"
#include "resources/library_manager.h"

LibraryManagerUnittest::LibraryManagerUnittest() {
}

int
LibraryManagerUnittest::test_get_motif() {
   
    MotifOP m  = LibraryManager::getInstance().get_motif("HELIX.IDEAL");
    MotifOP m1 = LibraryManager::getInstance().get_motif("TWOWAY.1GID.0");
    try {
        MotifOP m2 = LibraryManager::getInstance().get_motif("TWOWAY.1GID.1110");
        std::cout << "did not catch exception" << std::endl;
        exit(0);
    } catch(String const & e) {}

    MotifOP bp = LibraryManager::getInstance().get_motif("GC=GC");
    
    
    return 1;
}

int
LibraryManagerUnittest::run() {
    
    if (test_get_motif() == 0)             { std::cout << "test_get_motif failed" << std::endl; }

    return 0;
}