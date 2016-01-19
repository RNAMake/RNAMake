//
//  motif_state_tree_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_state_tree_unittest.h"
#include "build/build_motif_tree.h"

#include "motif_data_structures/motif_state_tree.h"
#include "resources/resource_manager.h"

namespace unittests {
namespace motif_structures {

int
MotifStateTreeUnittest::test_creation() {
    MotifStateTree mst;
    auto sterics = mst.get_bool_option("sterics");
    if(sterics != true) {
        throw UnittestException("was unable to retreieve correct option value");
    }
    
    mst.set_option_value("sterics", false);
    sterics = mst.get_bool_option("sterics");
    if(sterics != false) {
        throw UnittestException("was unable to retreieve correct option value");
    }
   return 1;
}

int
MotifStateTreeUnittest::test_add_state() {
    auto ms = ResourceManager::getInstance().get_state("HELIX.IDEAL");
    MotifStateTree mst;
    mst.add_state(ms);
    mst.add_state(ms);
    if(mst.size() != 2) { return 0; }
    return 1;
}

int
MotifStateTreeUnittest::test_from_mt() {
    BuildMotifTree builder;
    
    for(int i = 0; i < 50; i++) {
    
        auto mt = builder.build(10);
        MotifStateTree mst;
        mst.setup_from_mt(mt);
        if(mt->size() != mst.size()) { return 0;}
        
        auto d1 = mt->last_node()->data()->ends()[0]->d();
        auto d2 =  mst.last_node()->data()->cur_state->end_states()[0]->d();
        
        if(d1.distance(d2) > 4) { return 0; }
    
    }
    
    return 1;
}

int
MotifStateTreeUnittest::test_to_motif_tree() {
    BuildMotifTree builder;
    auto mt = builder.build(10);
    MotifStateTree mst;
    mst.setup_from_mt(mt);

    auto mt2 = mst.to_motif_tree();
    
    return 1;
}

int
MotifStateTreeUnittest::test_replace_state() {
    BuildMotifTree builder;
    auto mt = builder.build(10);
    MotifStateTree mst;
    mst.setup_from_mt(mt);
    auto state = ResourceManager::getInstance().get_state("HELIX.IDEAL.20");
    mst.replace_state(4, state);
    return 1;

}

int
MotifStateTreeUnittest::run() {
    if (test_creation() == 0)      { std::cout << "test_creation failed" << std::endl;  }
    if (test_add_state() == 0)     { std::cout << "test_add_state failed" << std::endl;  }
    //if (test_from_mt() == 0)        { std::cout << "test_from_mt failed" << std::endl;  }
    //if (test_to_motif_tree() == 0)  { std::cout << "test_to_motif_tree failed" << std::endl;  }
    if (test_replace_state() == 0)      { std::cout << "test_replace_state failed" << std::endl; }
    
    return 1;
}

int
MotifStateTreeUnittest::run_all() {
    String name = "MotifStateTreeUnittest";
    typedef int (MotifStateTreeUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_creation"      ] = &MotifStateTreeUnittest::test_creation;
    func_map["test_add_state"     ] = &MotifStateTreeUnittest::test_add_state;
    func_map["test_from_mt"       ] = &MotifStateTreeUnittest::test_from_mt;
    func_map["test_to_motif_tree" ] = &MotifStateTreeUnittest::test_to_motif_tree;
    func_map["test_replace_state" ] = &MotifStateTreeUnittest::test_replace_state;

    
    int failed = 0;
    for(auto const & kv : func_map) {
        try {
            (this->*kv.second)();
        }
        catch(std::exception const & e) {
            std::cout << name << "::" << kv.first << " returned ERROR! : " << e.what() << std::endl;
            failed++;
        }
    }
    return failed;
}
    
}
}
