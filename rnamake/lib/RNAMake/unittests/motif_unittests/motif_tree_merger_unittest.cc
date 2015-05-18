//
//  motif_tree_merger_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_merger_unittest.h"
#include "motif/motif_tree.h"
#include "resources/motif_library.h"

int
MotifTreeMergerUnittest::test_merger() {
    MotifLibrary mlib(HELIX);
    MotifTree mt;
    mt.add_motif(mlib.get_motif("HELIX.IDEAL"));
    mt.add_motif(mlib.get_motif("HELIX.IDEAL"));
    for(int i = 0; i < 10; i++) {
        mt.add_motif(mlib.get_motif("HELIX.IDEAL"));
    }
    
    PoseOP pose = mt.to_pose();
    if(pose->chains().size() != 2) { return 0; }
    return 1;
}

int
MotifTreeMergerUnittest::run() {
    if (test_merger() == 0)            { std::cout << "test_merger failed" << std::endl;  }
    return 0;
}

void
MotifTreeMergerUnittest::run_all() {
    String name = "MotifTreeMergerUnittest";
    typedef int (MotifTreeMergerUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_merger"   ] = &MotifTreeMergerUnittest::test_merger;

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
