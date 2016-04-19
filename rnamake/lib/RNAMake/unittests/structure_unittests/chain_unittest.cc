//
//  chain_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "is_equal.hpp"
#include "chain_unittest.h"
#include "math/numerical.h"
#include "util/file_io.h"

namespace unittests {

ChainUnittest::ChainUnittest() {
    
    auto path = unittest_resource_dir() + "/chain/test_str_to_chain.dat";
    auto lines = get_lines_from_file(path);
    
    rts_ = ResidueTypeSet();
    c_ = std::make_shared<Chain>(lines[0], rts_);
}

int
ChainUnittest::test_to_str() {

    auto s = c_->to_str();
    auto c2 = std::make_shared<Chain>(s, rts_);
    
    failUnless(are_chains_equal(c_, c2, 0), "chains should be equal");
    
    return 1;
}

    
int
ChainUnittest::test_copy() {
    auto c_copy = std::make_shared<Chain>(*c_);
    failUnless(are_chains_equal(c_, c_copy, 1), "chains should be equal");
    return 1;
    
}
    
int
ChainUnittest::test_subchain() {
    
    try {
        ChainOP sc1 = c_->subchain(-1, 5);
        std::cout << "did not catch exception\n" << std::endl;
        exit(0);
    } catch (char const * e) {}
    
    try {
        ChainOP sc1 = c_->subchain(1, 1);
        std::cout << "did not catch exception\n" << std::endl;
        exit(0);
    } catch (char const * e) {}
    
    try {
        ChainOP sc1 = c_->subchain(1, 10000);
        std::cout << "did not catch exception\n" << std::endl;
        exit(0);
    } catch (char const * e) {}
    
    ChainOP sc = c_->subchain(1, 5);
    if(sc->residues().size() != 4) { return 0; }
    
    return 1;
}

int
ChainUnittest::test_to_pdb() {
    int i = 1;
    String s = c_->to_pdb_str(i);
    
    return 1;
}

int
ChainUnittest::run() {
    
    test_to_str();
    test_subchain();
    test_to_pdb();
    test_copy();
    return 1;
}

int
ChainUnittest::run_all() {
    String name = "ChainUnittest";
    typedef int (ChainUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_to_str"   ] = &ChainUnittest::test_to_str;
    func_map["test_subchain" ] = &ChainUnittest::test_subchain;
    func_map["test_to_pdb"   ] = &ChainUnittest::test_to_pdb;
    
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
