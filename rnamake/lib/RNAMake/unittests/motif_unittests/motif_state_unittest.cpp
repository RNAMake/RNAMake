//
//  motif_state_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_state_unittest.h"
#include "structure/basepair_state.h"
#include "resources/resource_manager.h"

int
MotifStateUnittest::test_creation() {
    auto m = ResourceManager::getInstance().get_motif("HELIX.IDEAL");
    auto ms = m->get_state();
    
    return 1;
}

int
MotifStateUnittest::test_copy() {
    auto m = ResourceManager::getInstance().get_motif("HELIX.IDEAL");
    auto ms = m->get_state();
    auto ms_copy = std::make_shared<MotifState>(ms->copy());
    
    return 1;
    
}

int
MotifStateUnittest::test_to_str() {
    auto m = ResourceManager::getInstance().get_motif("HELIX.IDEAL");
    auto ms = m->get_state();
    auto s = ms->to_str();
    auto ms1 = str_to_motif_state(s);
    
    /*for(int i = 0; i < ms->end_states().size(); i++) {
        if(!are_BasepairStates_equal(*ms->end_states()[i], *ms1.end_states()[i])) { return 0; }
    }*/
    
    return 1;
}

int
MotifStateUnittest::test_align() {
    String path = base_dir() + "/rnamake/unittests/resources/motifs/tetraloop_receptor_min";
    ResourceManager::getInstance().add_motif(path);
    auto m = ResourceManager::getInstance().get_motif("tetraloop_receptor_min", "", "A228-A246");
    auto bp_state = m->ends()[1]->state();
    auto test_state = ResourceManager::getInstance().get_motif("HELIX.IDEAL.3")->get_state();
    auto d1 = bp_state->d();
    align_motif_state(bp_state, test_state);
    auto d2 = test_state->end_states()[0]->d();
    auto dist = d1.distance(d2);
    if(dist > 0.5) {
        return 0;
    }
    
    
    return 1;
}

int
MotifStateUnittest::run() {
    if (test_creation() == 0)        { std::cout << "test_creation failed" << std::endl;  }
    if (test_copy() == 0)            { std::cout << "test_copy failed" << std::endl;  }
    if (test_to_str() == 0)          { std::cout << "test_to_str failed" << std::endl;  }
    if (test_align() == 0)           { std::cout << "test_align failed" << std::endl;  }
    return 1;
}