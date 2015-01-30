//
//  main.cpp
//  basepair.unittest
//
//  Created by Joseph Yesselman on 1/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "structure.h"
#include "basepair.h"
#include "numerical.h"


Structure
get_test_structure() {
    String file = "test_str_to_structure.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    getline(input, line);
    Structure s = str_to_structure(line, rts);
    return s;
}

Strings
get_lines_from_file(String const fname) {
    String line;
    Strings lines;
    std::ifstream input;
    input.open(fname);
    while ( input.good() ) {
        getline(input, line);
        if( line.length() < 10 ) { break; }
        lines.push_back(line);

    }
    return lines;

}

int
test_creation() {
    Structure s = get_test_structure();
    Residue res1 = s.get_residue(103, "A", "");
    Residue res2 = s.get_residue(104, "A", "");
    Matrix r(0.0);
    Basepair bp (res1, res2, r, "c...");
    BasepairState bpstate = bp.state();
    return 1;
}

int
test_str_to_basepair_state() {
    String file = "test_str_to_basepairstate.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    while ( input.good() ) {
        getline(input, line);
        if( line.length() < 10 ) { break; }
        BasepairState bpstate = str_to_basepairstate(line);
    }
    return 1;
}

int
test_basepair_state_to_str() {
    String file = "test_str_to_basepairstate.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    while ( input.good() ) {
        getline(input, line);
        if( line.length() < 10 ) { break; }
        BasepairState bpstate = str_to_basepairstate(line);
        String str = bpstate.to_str();
        BasepairState bpstate2 = str_to_basepairstate(line);
        if (! (are_BasepairStates_equal(bpstate, bpstate2) )) { return 0; }

    }
    return 1;
}

int
test_get_transforming_r_and_t_test() {
    Strings lines = get_lines_from_file("get_transforming_r_and_t_test.dat");
    BasepairState dummy;
    for (auto const & l : lines) {
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


int main(int argc, const char * argv[]) {
    if (test_creation() == 0)                            { std::cout << "test_creation failed" << std::endl; }
    if (test_str_to_basepair_state() == 0)               { std::cout << "test_str_to_basepair_state failed" << std::endl; }
    if (test_basepair_state_to_str() == 0)               { std::cout << "test_basepair_state_to_str failed" << std::endl; }
    if (test_get_transforming_r_and_t_test() == 0)       { std::cout << "test_get_transforming_r_and_t_test failed" << std::endl; }

    return 0;
}









