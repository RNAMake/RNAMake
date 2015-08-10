//
//  added_motif_library_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/9/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "added_motif_library_unittest.h"
#include "resources/resource_manager.h"
#include "resources/added_motif_library.h"

int
AddedMotifLibraryUnittest::test_add_motif() {
    AddedMotifLibrary mlib;
    mlib.add_motif(ResourceManager::getInstance().get_motif("HELIX.IDEAL"));
    
    return 1;
}

int
AddedMotifLibraryUnittest::test_get() {
    AddedMotifLibrary mlib;
    mlib.add_motif(ResourceManager::getInstance().get_motif("HELIX.IDEAL"));
    
    auto m = mlib.get("HELIX.IDEAL");
    auto m1 = mlib.get("", "CC_LL_GG_RR");
    auto m2 = mlib.get("", "", "A5-B7");
    auto m3 = mlib.get("HELIX.IDEAL", "CC_LL_GG_RR");
    auto m4 = mlib.get("HELIX.IDEAL", "", "A5-B7");
    auto m5 = mlib.get("HELIX.IDEAL", "CC_LL_GG_RR", "A5-B7");
    return 1;
}

int
AddedMotifLibraryUnittest::test_get_multi() {
    AddedMotifLibrary mlib;
    mlib.add_motif(ResourceManager::getInstance().get_motif("HELIX.IDEAL"));
    mlib.add_motif(ResourceManager::getInstance().get_motif("HELIX.IDEAL"));
    auto motifs = mlib.get_multi("HELIX.IDEAL");
    if(motifs.size() != 2) { return 0; }
    return 1;
}

int
AddedMotifLibraryUnittest::test_contains() {
    AddedMotifLibrary mlib;
    mlib.add_motif(ResourceManager::getInstance().get_motif("HELIX.IDEAL"));
    if(!mlib.contains("HELIX.IDEAL")) { return 0; }
    if(mlib.contains("TEST"))         { return 0; }
    if(mlib.contains("", "TEST"))     { return 0; }
    return 1;
}

int
AddedMotifLibraryUnittest::run() {
    if (test_add_motif() == 0)        { std::cout << "test_add_motif failed" << std::endl; }
    if (test_get() == 0)              { std::cout << "test_get failed" << std::endl; }
    if (test_get_multi() == 0)        { std::cout << "test_get failed" << std::endl; }
    if (test_contains() == 0)         { std::cout << "test_contain failed" << std::endl; }
    return 0;
}