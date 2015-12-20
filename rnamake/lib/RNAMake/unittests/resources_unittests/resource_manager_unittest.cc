//
//  resource_manager_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/9/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "resource_manager_unittest.h"

#include "util/settings.h"
#include "resources/resource_manager.h"

int
ResourceManagerUnittest::test_get_motif() {
    auto m = ResourceManager::getInstance().get_motif("HELIX.IDEAL");
    if(m->name() != "HELIX.IDEAL") { return 0; }
    return 1;
}

int
ResourceManagerUnittest::test_add_motif() {
    auto path = base_dir() + "/rnamake/unittests/resources/motifs/fmn_min";
    ResourceManager::getInstance().add_motif(path);
    
    auto m = ResourceManager::getInstance().get_motif("fmn_min", "", "A15-A22");

    //wrong end_name
    try {
        auto m1 =  ResourceManager::getInstance().get_motif("fmn_min", "", "A16-A22");
        throw std::runtime_error("got unexpected error");
    }
    catch(ResourceManagerException e) {}
    catch(...) { return 0; }
    
    return 1;
}

void
ResourceManagerUnittest::test_memory() {
    auto names = Strings{"HELIX.IDEAL", "HELIX.IDEAL.2", "HELIX.IDEAL.3", "HELIX.IDEAL.4", "HELIX.IDEAL.5"};
    auto rng = RandomNumberGenerator();
    
    auto mlib = std::make_shared<MotifSqliteLibrary>("ideal_helices");
    int count = 0;
    for(int i = 0; i < 1000000; i++) {
        //auto m = mlib->get(names[rng.randrange(names.size())]);
        auto m = ResourceManager::getInstance().get_motif(names[rng.randrange(names.size())]);
        count += m->residues().size();
    }
}


int
ResourceManagerUnittest::run() {
    //if (test_get_motif() == 0)        { std::cout << "test_get_motif failed" << std::endl; }
    //if (test_add_motif() == 0)        { std::cout << "test_add_motif failed" << std::endl; }
    test_memory();
    
    return 0;
}