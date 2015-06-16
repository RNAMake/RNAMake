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
LibraryManagerUnittest::test_add_motif() {
    String m_path = base_dir() + "/rnamake/unittests/resources/motifs/fmn_min";
    LibraryManager::getInstance().add_motif(m_path);
    
    MotifOP m  = LibraryManager::getInstance().get_motif("fmn_min");
    
    return 1;
}

int
LibraryManagerUnittest::run() {
    
    if (test_get_motif() == 0)             { std::cout << "test_get_motif failed" << std::endl; }
    if (test_add_motif() == 0)             { std::cout << "test_add_motif failed" << std::endl; }
    
    return 0;
}

void
LibraryManagerUnittest::run_all() {
    String name = "LibraryManagerUnittest";
    typedef int (LibraryManagerUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_get_motif" ] = &LibraryManagerUnittest::test_get_motif;
    func_map["test_add_motif" ] = &LibraryManagerUnittest::test_add_motif;

    for(auto const & kv : func_map) {
        try {
            int result = (this->*kv.second)();
            if(result == 0) {
                std::cout << name << "::" << kv.first << " FAILED!" << std::endl;
            }
        }
        catch(...) {
            std::cout << name << "::" << kv.first << " returned ERROR!" << std::endl;
        }
        
    }
}