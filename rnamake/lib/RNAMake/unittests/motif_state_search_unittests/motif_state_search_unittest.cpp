//
//  motif_state_search_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_state_search_unittest.h"
#include "build/build_motif_tree.h"



namespace unittests {
namespace motif_state_search {

void
MotifStateSearchUnittest::test_creation() {
    MotifStateSearch mss;
}

void
MotifStateSearchUnittest::test_search() {
    MotifStateSearch mss;
    mss.set_option_value("accept_score", 10);
    mss.set_option_value("max_node_level", 20);
    mss.set_option_value("sterics", false);
    BuildMotifTree builder;
    auto mt = builder.build(10);
    auto start = mt->get_node(0)->data()->ends()[0]->state();
    auto end   = mt->last_node()->data()->ends()[1]->state();
    
    //mt->to_pdb("test.pdb");
    
    mss.setup(start, end);
    auto sol = mss.next();

    if(sol == nullptr) {
        throw UnittestException("could not find a solution to search");
    }
    mt->to_pdb("test.pdb", 1);
    sol->to_pdb("solution.pdb", 1);
}

/*int
MotifStateSearchUnittest::test_aligner() {
    BuildMotifTree builder;
    auto mt = builder.build(10);
    RandomNumberGenerator rng;
    
    MotifStateOP cur = mt->get_node(0)->data()->get_state();
    MotifStateOPs states;
    for(auto const & n : *mt) {
        states.push_back(n->data()->get_state());
    }
    aligner a;
    auto start = mt->get_node(0)->data()->ends()[0]->state();
    MotifStateOP ref_state;
    
    for(int i = 0; i < 10000000; i++) {
        ref_state   = states[rng.randrange(9)];
        a.get_aligned_motif_state(start, cur, ref_state);
    }
    
    
    
    return 1;
}*/

int
MotifStateSearchUnittest::run() {
    test_creation();
    test_search();
    //if(test_creation() == 0)       { std::cout << "test_creation failed" << std::endl; }
    //if(test_search() == 0)         { std::cout << "test_search failed" << std::endl; }
    //if(test_aligner() == 0)        { std::cout << "test_aligner failed" << std::endl; }
    return 1;
}
    

int
MotifStateSearchUnittest::run_all() {
    String name = "MotifStateSearchUnittest";
    typedef void (MotifStateSearchUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_creation"   ] = &MotifStateSearchUnittest::test_creation;
    func_map["test_search"     ] = &MotifStateSearchUnittest::test_search;
    
    
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

