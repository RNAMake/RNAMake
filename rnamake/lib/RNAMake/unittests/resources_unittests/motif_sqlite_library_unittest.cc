//
//  motif_sqlite_library_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_sqlite_library_unittest.h"

int
MotifSqliteLibraryUnittest::test_creation() {
    MotifSqliteLibrary mlib("bp_steps");
    return 1;
}

int
MotifSqliteLibraryUnittest::test_get() {
    MotifSqliteLibrary mlib("ideal_helices");
    auto m = mlib.get("HELIX.IDEAL");
    std::cout << m->name() << std::endl;
    return 1;
}

int
MotifSqliteLibraryUnittest::run() {
    if (test_creation() == 0)        { std::cout << "test_creations failed" << std::endl; }
    if (test_get() == 0)             { std::cout << "test_get failed" << std::endl; }

    return 0;
    
}