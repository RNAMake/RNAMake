//
//  motif_tree_state_library_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/15/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_state_library_unittest.h"
#include "util/settings.h"
#include "motif_tree_state/motif_tree_state_library.h"

int
MotifTreeStateLibraryUnittest::test_creation() {
    String m_path = base_dir() + "/rnamake/unittests/test_twoway.new.me";
    MotifTreeStateLibrary mts_lib(m_path, 1);
    return 1;
}


int
MotifTreeStateLibraryUnittest::run() {
    if (test_creation() == 0)          { std::cout << "test_creation failed" << std::endl;  }
    return 1;
}

void
MotifTreeStateLibraryUnittest::run_all() {
    String name = "MotifTreeStateLibraryUnittest";
    typedef int (MotifTreeStateLibraryUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_creation"   ] = &MotifTreeStateLibraryUnittest::test_creation;
    
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
