//
//  motif_state_ensemble_sqlite_library_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_state_ensemble_sqlite_library_unittest.h"

#include "motif_state_sqlite_library_unittest.h"
#include "resources/motif_state_ensemble_sqlite_library.h"

int
MotifStateEnsembleSqliteLibraryUnittest::test_creation() {
    MotifStateEnsembleSqliteLibrary mse_lib("bp_steps");
    return 1;
}

int
MotifStateEnsembleSqliteLibraryUnittest::test_get() {
    MotifStateEnsembleSqliteLibrary mse_lib("bp_steps");
    auto m  = mse_lib.get("GG_LL_CC_RR");
    auto m1 = mse_lib.get("", "1");
    return 1;
}

int
MotifStateEnsembleSqliteLibraryUnittest::test_get_random() {
    MotifStateEnsembleSqliteLibrary mlib("bp_steps");
    auto m = mlib.get_random();
    for(int i = 0; i < 100; i++) {
        m = mlib.get_random();
    }
    
    return 1;
}

int
MotifStateEnsembleSqliteLibraryUnittest::test_all() {
    MotifStateEnsembleSqliteLibrary mlib("bp_steps");
    mlib.load_all();
    
    int count = 0;
    for(auto const & m : mlib) {
        auto m1 = m;
        count ++;
    }
    if(count == 0) { return 0; }
    return 1;
}

int
MotifStateEnsembleSqliteLibraryUnittest::test_contains() {
    MotifStateEnsembleSqliteLibrary mlib("bp_steps");
    if(!mlib.contains("GG_LL_CC_RR")) { return 0; }
    if(mlib.contains("TEST"))         { return 0; }
    if(mlib.contains("", "TEST"))     { return 0; }
    
    return 1;
}


int
MotifStateEnsembleSqliteLibraryUnittest::run() {
    if (test_creation() == 0)        { std::cout << "test_creations failed" << std::endl; }
    if (test_get() == 0)             { std::cout << "test_get failed" << std::endl; }
    if (test_get_random() == 0)      { std::cout << "test_get_random failed" << std::endl; }
    if (test_all() == 0)             { std::cout << "test_all failed" << std::endl; }
    if (test_contains() == 0)        { std::cout << "test_contains failed" << std::endl; }
    return 0;
    
}