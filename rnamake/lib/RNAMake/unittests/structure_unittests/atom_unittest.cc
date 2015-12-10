//
//  atom_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <functional>
#include <map>

#include "base/types.h"
#include "atom_unittest.h"


namespace unittests {
    
int
AtomUnittest::test_creation() {
    Atom atom("P", Point(0, 1, 2));
    Point p(0, 1, 2);
    if (!are_xyzVector_equal(atom.coords(), p)) {
        return 0;
    }
    return 1;
}

int
AtomUnittest::test_to_pdb_str() {
    Atom atom("H1", Point(1, 2, 3));
    String s = atom.to_pdb_str(1);
    String ref("ATOM      1  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P\n");
    if (s.compare(ref) != 0) {
        return 1;
    }
    String s2 = atom.to_pdb_str(10);
    String ref2("ATOM     10  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P\n");
    if (s2.compare(ref2) != 0) {
        return 1;
    }
    
    return 1;
}

int
AtomUnittest::test_str_to_atom() {
    auto a = Atom("H1", Point(1, 2, 3));
    String s = a.to_str();
    auto a2 = Atom(s);
    if(not are_xyzVector_equal(a.coords(), a2.coords())) {
        return 0;
    }
    if(not (a.name().compare(a2.name()) == 0)) {
        return 0;
    }
    return 1;
}

int
AtomUnittest::test_copy() {
    auto a  = Atom("H1", Point(1, 2, 3));
    auto a2 = Atom(a);
    if(not are_xyzVector_equal(a.coords(), a2.coords())) {
        return 0;
    }
    if(not (a.name().compare(a2.name()) == 0)) {
        return 0;
    }
    return 1;
}

int
AtomUnittest::run() {
    
    if (test_creation() == 0)    {  std::cout << "test_creation failed" << std::endl; }
    if (test_to_pdb_str() == 0)  {  std::cout << "test_to_pdb_str failed" << std::endl; }
    if (test_str_to_atom() == 0) {  std::cout << "test_str_to_atom failed" << std::endl; }
    if (test_copy() == 0)        {  std::cout << "test_copy failed" << std::endl; }
    
    return 1;
}

int
AtomUnittest::run_all() {
    String name = "AtomUnittest";
    typedef int (AtomUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_creation"   ] = &AtomUnittest::test_creation;
    func_map["test_to_pdb_str" ] = &AtomUnittest::test_to_pdb_str;
    func_map["test_str_to_atom"] = &AtomUnittest::test_str_to_atom;
    func_map["test_copy"       ] = &AtomUnittest::test_copy;
    
    int failed = 0;
    for(auto const & kv : func_map) {
        try {
            (this->*kv.second)();
        }
        catch(...) {
            std::cout << name << "::" << kv.first << " returned ERROR!" << std::endl;
            failed += 1;
        }
        
    }
    
    return failed;
}


}








