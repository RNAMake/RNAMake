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
MotifSqliteLibraryUnittest::test_get_random() {
    MotifSqliteLibrary mlib("ideal_helices");
    auto m = mlib.get_random();
    for(int i = 0; i < 100; i++) {
        m = mlib.get_random();
    }
    
    return 1;
}

int
MotifSqliteLibraryUnittest::test_get_multi() {
    MotifSqliteLibrary mlib("twoway");
    auto m = mlib.get_random();
    auto motifs_1 = mlib.get_multi(m->name());
    auto motifs_2 = mlib.get_multi("", m->end_ids()[0]);
    return 1;
}

int
MotifSqliteLibraryUnittest::test_all() {
    MotifSqliteLibrary mlib("ideal_helices");
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
MotifSqliteLibraryUnittest::test_contains() {
    MotifSqliteLibrary mlib("ideal_helices");
    if(!mlib.contains("HELIX.IDEAL")) { return 0; }
    if(mlib.contains("TEST"))         { return 0; }
    if(mlib.contains("", "TEST"))     { return 0; }
    
    return 1;
}

void
MotifSqliteLibraryUnittest::test_memory() {
    MotifSqliteLibrary mlib("ideal_helices");
    auto names = Strings{"HELIX.IDEAL", "HELIX.IDEAL.2", "HELIX.IDEAL.3", "HELIX.IDEAL.4"};
    auto rng = RandomNumberGenerator();
    int count = 0;
    for(int i = 0; i < 1000000; i++) {
        int in = mlib.contains(names[rng.randrange(names.size())]);
        count += in;
    }
    
}


int
MotifSqliteLibraryUnittest::run() {
    /*if (test_creation() == 0)        { std::cout << "test_creations failed" << std::endl; }
    if (test_get() == 0)             { std::cout << "test_get failed" << std::endl; }
    if (test_get_random() == 0)      { std::cout << "test_get_random failed" << std::endl; }
    if (test_all() == 0)             { std::cout << "test_all failed" << std::endl; }
    if (test_get_multi() == 0)       { std::cout << "test_get_multi failed" << std::endl; }
    if (test_contains() == 0)        { std::cout << "test_contains failed" << std::endl; }*/
    test_memory();
    return 0;
    
}