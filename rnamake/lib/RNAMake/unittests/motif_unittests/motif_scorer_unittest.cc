//
//  motif_scorer_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_scorer_unittest.h"
#include "motif/motif_scorer.h"
#include "resources/resource_manager.h"

int
MotifScorerUnittest::test_score() {
    MotifScorer ms;
    MotifOP m  = ResourceManager::getInstance().get_motif("HELIX.IDEAL");
    float score = ms.score(m);
    if(score - -3.4 > 0.1) { return 1; }
    
    return 1;
}



int
MotifScorerUnittest::run() {
    if (test_score() == 0)            { std::cout << "test_score failed" << std::endl;  }
    return 0;
}

void
MotifScorerUnittest::run_all() {
    String name = "MotifScorerUnittest";
    typedef int (MotifScorerUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_score"   ] = &MotifScorerUnittest::test_score;
    
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
