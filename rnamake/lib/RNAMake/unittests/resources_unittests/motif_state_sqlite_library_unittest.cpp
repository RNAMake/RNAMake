//
//  motif_state_sqlite_library.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_state_sqlite_library_unittest.h"
#include "resources/motif_state_sqlite_library.h"

int
MotifStateSqliteLibraryUnittest::test_creation() {
    MotifStateSqliteLibrary mlib("bp_steps");
    return 1;
}

int
MotifStateSqliteLibraryUnittest::test_get() {
    MotifStateSqliteLibrary mlib("ideal_helices");
    auto m  = mlib.get("HELIX.IDEAL");
    auto m1 = mlib.get("", "CC_LL_GG_RR");
    auto m2 = mlib.get("", "", "A5-B7");
    auto m3 = mlib.get("HELIX.IDEAL", "CC_LL_GG_RR");
    auto m4 = mlib.get("HELIX.IDEAL", "", "A5-B7");
    auto m5 = mlib.get("HELIX.IDEAL", "CC_LL_GG_RR", "A5-B7");
    auto m6 = mlib.get("", "", "", "1");
    return 1;
}

int
MotifStateSqliteLibraryUnittest::test_get_random() {
    MotifStateSqliteLibrary mlib("ideal_helices");
    auto m = mlib.get_random();
    for(int i = 0; i < 100; i++) {
        m = mlib.get_random();
    }
    
    return 1;
}

int
MotifStateSqliteLibraryUnittest::test_get_multi() {
    MotifStateSqliteLibrary mlib("twoway");
    auto m = mlib.get_random();
    auto motifs_1 = mlib.get_multi(m->name());
    auto motifs_2 = mlib.get_multi("", m->end_ids()[0]);
    return 1;
}

int
MotifStateSqliteLibraryUnittest::test_all() {
    MotifStateSqliteLibrary mlib("ideal_helices");
    mlib.load_all();
    
    int count = 0;
    for(auto const & m : mlib) {
        auto m1 = m;
        count ++;
    }
    if(count != 21) { return 0;}
    return 1;
}

int
MotifStateSqliteLibraryUnittest::test_contains() {
    MotifStateSqliteLibrary mlib("ideal_helices");
    if(!mlib.contains("HELIX.IDEAL")) { return 0; }
    if(mlib.contains("TEST"))         { return 0; }
    if(mlib.contains("", "TEST"))     { return 0; }
    
    return 1;
}


int
MotifStateSqliteLibraryUnittest::run() {
    if (test_creation() == 0)        { std::cout << "test_creations failed" << std::endl; }
    if (test_get() == 0)             { std::cout << "test_get failed" << std::endl; }
    if (test_get_random() == 0)      { std::cout << "test_get_random failed" << std::endl; }
    if (test_all() == 0)             { std::cout << "test_all failed" << std::endl; }
    if (test_get_multi() == 0)       { std::cout << "test_get_multi failed" << std::endl; }
    if (test_contains() == 0)        { std::cout << "test_contains failed" << std::endl; }
    return 0;
    
}