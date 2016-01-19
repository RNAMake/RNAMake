//
//  chain_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "chain_unittest.h"
#include "math/numerical.h"
#include "util/file_io.h"

namespace unittests {

ChainUnittest::ChainUnittest() {
    
    String path = unittest_resource_dir() + "/chain/test_str_to_chain.dat";
    Strings lines = get_lines_from_file(path);
    
    ResidueTypeSet rts;
    c_ = std::make_shared<Chain>(lines[0], rts);
}

int
ChainUnittest::test_to_str() {

    String s = c_->to_str();
    ResidueTypeSet rts;
    auto c2 = Chain(s, rts);
    
    ResidueOPs org_res = c_->residues();
    ResidueOPs new_res = c2.residues();
    
    for(int i = 0; i < org_res.size(); i++) {
        if(org_res[i]->name().compare(new_res[i]->name()) != 0) { return 0; }
        if(org_res[i]->chain_id().compare(new_res[i]->chain_id())) { return  0; }
        if(org_res[i]->num() != new_res[i]->num()) { return 0; }
        
        AtomOPs atoms_1 = org_res[i]->atoms();
        AtomOPs atoms_2 = new_res[i]->atoms();
        
        for(int j = 0; j < atoms_1.size(); j++) {
            if(atoms_1[j] == nullptr && atoms_2[j] != nullptr) { return 0; }
            if(atoms_1[j] != nullptr && atoms_2[j] == nullptr) { return 0; }
            if(atoms_1[j] == nullptr && atoms_2[j] == nullptr) { continue; }
            
            if(!are_xyzVector_equal(atoms_1[j]->coords(), atoms_2[j]->coords())) { return 0; }
        }
    }
    
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
    
    if (test_to_str() == 0)             { std::cout << "test_to_str failed" << std::endl; }
    if (test_subchain() == 0)           { std::cout << "test_subchain failed" << std::endl; }
    if (test_to_pdb() == 0)             { std::cout << "test_to_pdb failed" << std::endl; }
    
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
