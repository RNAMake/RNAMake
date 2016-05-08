//
//  atom_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <functional>
#include <map>

#include "is_equal.hpp"
#include "base/types.h"
#include "atom_unittest.h"


namespace unittests {
    
int
AtomUnittest::test_creation() {
    auto a = std::make_shared<Atom>("P", Point(0, 1, 2));
    auto p = Point(0, 1, 2);
 
    failUnless(a->coords() == p, "did properly create atomic coordinates");
    
    return 1;
}

int
AtomUnittest::test_to_pdb_str() {
    auto atom = std::make_shared<Atom>("P", Point(1, 2, 3));
    auto s = atom->to_pdb_str(1);
    auto ref = String("ATOM      1  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P\n");
    
    failUnless(s == ref, "did not produce right string output");
    
    auto s2 = atom->to_pdb_str(10);
    auto ref2 = String("ATOM     10  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P\n");
    
    failUnless(s2 == ref2, "did not produce right string output");

    return 1;
}

int
AtomUnittest::test_str_to_atom() {
    auto a = std::make_shared<Atom>("P", Point(0, 1, 2));
    auto s = a->to_str();
    auto a2 = std::make_shared<Atom>(s);
    
    failUnless(are_atoms_equal(a, a2), "atoms are not equal");
    
    return 1;
}

int
AtomUnittest::test_copy() {
    auto a = std::make_shared<Atom>("P", Point(0, 1, 2));
    auto a2 = std::make_shared<Atom>(*a);
   
    failUnless(are_atoms_equal(a, a2), "atoms are not equal but should be");
    
    return 1;
}
    
int
AtomUnittest::test_to_str() {
    auto a = std::make_shared<Atom>("H1", Point(0, 1, 2));
    auto s = a->to_str();
    failUnless(s != "H1 0.0 1.0 2.0", "did not get correct string");
    
    return 1;
}

int
AtomUnittest::run() {
    test_creation();
    test_to_pdb_str();
    test_str_to_atom();
    test_copy();
    test_to_str();

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
    func_map["test_to_str"     ] = &AtomUnittest::test_to_str;

    
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








