//
//  motif_tree_merger_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_merger_unittest.h"
#include "motif/motif_tree.h"
#include "resources/resource_manager.h"

int
MotifTreeMergerUnittest::test_merger() {
    MotifOP m = ResourceManager::getInstance().get_motif("HELIX.IDEAL");
    MotifTree mt;
    for(int i = 0; i < 2; i++) { mt.add_motif(m);}
    
    PoseOP p = mt.to_pose();
    if(p->chains().size() != 2) { return 0; }
    
    auto motifs = p->motifs(MotifType::ALL);
    auto ss_motifs = p->secondary_structure()->motifs("ALL");
    
    if(motifs.size() != ss_motifs.size()) { return 0; }
 
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
