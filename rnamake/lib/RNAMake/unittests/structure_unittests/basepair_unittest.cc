//
//  basepair_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "basepair_unittest.h"
#include "util/file_io.h"
#include "structure/basepair_state.h"

BasepairUnittest::BasepairUnittest() {
    
    String path = unittest_resource_dir() + "/structure/test_str_to_structure.dat";
    Strings lines = get_lines_from_file(path);
    
    ResidueTypeSet rts;
    s_ = Structure(lines[0], rts);
    
}

int
BasepairUnittest::test_creation() {
    ResidueOP res1 = s_.get_residue(103, "A", "");
    ResidueOP res2 = s_.get_residue(104, "A", "");
    Matrix r(0.0);
    Basepair bp (res1, res2, r, "c...");
    auto bpstate = bp.state();
    return 1;
}

int
BasepairUnittest::test_str_to_basepair_state() {
    String path = unittest_resource_dir() + "/basepair/test_str_to_basepairstate.dat";
    Strings lines = get_lines_from_file(path);

    for(auto const & l : lines) {
        if(l.length() < 5) { break;}
        BasepairState bpstate = str_to_basepairstate(l);
    }
    return 1;
}

int
BasepairUnittest::test_basepair_state_to_str() {
    
    String path = unittest_resource_dir() + "/basepair/test_str_to_basepairstate.dat";
    Strings lines = get_lines_from_file(path);
    
    for(auto const & l : lines) {
        if(l.length() < 5) { break;}
        BasepairState bpstate = str_to_basepairstate(l);
        String str = bpstate.to_str();
        BasepairState bpstate2 = str_to_basepairstate(str);
        //if (! (are_BasepairStates_equal(bpstate, bpstate2) )) { return 0; }
        
    }
    return 1;
}

int
BasepairUnittest::test_get_transforming_r_and_t_test() {
    String path = unittest_resource_dir() + "/basepair/get_transforming_r_and_t_test.dat";
    Strings lines = get_lines_from_file(path);

    BasepairState dummy;
    for (auto const & l : lines) {
        if(l.length() < 5) { break;}
        Strings spl = split_str_by_delimiter(l, "|");
        BasepairState bpstate1 = str_to_basepairstate(spl[0]);
        BasepairState bpstate2 = str_to_basepairstate(spl[1]);
        Point t = vector_from_str(spl[2]);
        bpstate1.get_transforming_r_and_t(bpstate2, dummy);
        Point result = t + bpstate1.d();
        if (!(are_xyzVector_equal(result, dummy.d()))) {
            return 0;
        }
    }
    return 1;
}

int
BasepairUnittest::test_move() {
    ResidueOP res1 = s_.get_residue(103, "A", "");
    ResidueOP res2 = s_.get_residue(104, "A", "");
    Matrix r(0.0);
    Basepair bp (res1, res2, r, "c...");
    Point org_center = bp.d();
    s_.move(Point(10,10,10));
    Point new_center = bp.d();
    float dist = org_center.distance(new_center);
    if(dist - 30 > 0.001) { return 0; }
    return 1;
}

int
BasepairUnittest::run() {
    
    if (test_creation() == 0)                            { std::cout << "test_creation failed" << std::endl; }
    if (test_str_to_basepair_state() == 0)               { std::cout << "test_str_to_basepair_state failed" << std::endl; }
    /*if (test_basepair_state_to_str() == 0)               { std::cout << "test_basepair_state_to_str failed" << std::endl; }
    if (test_get_transforming_r_and_t_test() == 0)       { std::cout << "test_get_transforming_r_and_t_test failed" << std::endl; }
    if (test_move() == 0)                                { std::cout << "test_move failed" << std::endl; }
    */return 0;

    
}

int
BasepairUnittest::run_all() {
    String name = "BasepairUnittest";
    typedef int (BasepairUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_creation"                     ] = &BasepairUnittest::test_creation;
    func_map["test_str_to_basepair_state"        ] = &BasepairUnittest::test_str_to_basepair_state;
    func_map["test_get_transforming_r_and_t_test"] = &BasepairUnittest::test_get_transforming_r_and_t_test;
    func_map["test_move"                         ] = &BasepairUnittest::test_move;
    
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
    
    return 0;
}
