//
//  library_manager_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "library_manager_unittest.h"


LibraryManagerUnittest::LibraryManagerUnittest() {
    lm_ = LibraryManager();
    
}

int
LibraryManagerUnittest::test_get_motif() {
   
    MotifOP m  = lm_.get_motif("HELIX.IDEAL");
    MotifOP m1 = lm_.get_motif("TWOWAY.1GID.0");
    
    try {
        MotifOP m2 = lm_.get_motif("TWOWAY.1GID.1110");
        std::cout << "did not catch exception" << std::endl;
        exit(0);
    } catch(String const & e) {}

    
    return 1;
}

int
LibraryManagerUnittest::run() {
    
    if (test_get_motif() == 0)             { std::cout << "test_get_motif failed" << std::endl; }

    return 0;
}