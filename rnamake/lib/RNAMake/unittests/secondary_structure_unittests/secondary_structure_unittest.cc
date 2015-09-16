//
//  secondary_structure_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <memory.h>

#include "secondary_structure_unittest.h"

#include "util/uuid.h"
#include "secondary_structure/residue.h"
#include "secondary_structure/chain.h"
#include "secondary_structure/secondary_structure.h"
#include "secondary_structure/secondary_structure_factory.h"


int
SecondaryStructureUnittest::test_creation() {
    sstruct::Chain c;
    try {
        c.first();
        throw std::runtime_error("unexpected error");
    }
    catch(sstruct::SecondaryStructureException e) { }
    catch(...) { return 0; }
    
    sstruct::ResidueOPs res;
    res.push_back(std::make_shared<sstruct::Residue>("G", ".", 1, "A", Uuid()));
    res.push_back(std::make_shared<sstruct::Residue>("A", "(", 2, "A", Uuid()));
    res.push_back(std::make_shared<sstruct::Residue>("U", ".", 3, "A", Uuid()));
    sstruct::Chain c1(res);
    
    if(c1.sequence() != "GAU") {
        return 0;
    }
    
    String s = c1.to_str();
    auto c2 = sstruct::str_to_chain(s);

    if(c2.dot_bracket() != ".(.") {
        return 0;
    }
    
    auto c3 = c2.copy();
    if(c3.sequence() != "GAU") {
        return 0;
    }
    
    sstruct::SecondaryStructure ss("GGGG&CCCC", "((((&))))");
    if(ss.sequence() != "GGGG&CCCC") {
        return 0;
    }
    
    
    return 1;
}

int
SecondaryStructureUnittest::test_creation_residue() {
    sstruct::Residue r("G", ".", 1, "A", Uuid());
    String s = r.to_str();
    auto r2 = sstruct::str_to_residue(s);
    
    if(r.num() != r2.num() || r.chain_id() != r2.chain_id()) {
        return 0;
    }
    
    auto r3 = r.copy();
    
    if(!(r.uuid() == r3.uuid())) {
        return 0;
    }
    
    return 1;
}

int
SecondaryStructureUnittest::test_get_residue() {
    sstruct::SecondaryStructure ss("GGGG&CCCC", "((((&))))");
    auto r = ss.get_residue(0, "A");
    if(r == nullptr) { return 0; }
    if(r->num() != 0) { return 0; }
    
    auto r1 = ss.get_residue(1, "B");
    if(r1 != nullptr) { return 0; }
    
    auto r2 = ss.get_residue(14, "A");
    if(r2 != nullptr) { return 0; }

    auto r3 = ss.get_residue(r->uuid());
    if(r3 == nullptr) { return 0; }
        
    return 1;
}

int
SecondaryStructureUnittest::test_assign_end_id() {
    sstruct::SecondaryStructureFactory sf;
    auto ss = sf.get_structure("GGGG&CCCC", "((((&))))");
    auto end_id = sstruct::assign_end_id(ss, ss->ends()[0]);
    //std::cout << end_id << std::endl;
    auto end_id2 = sstruct::assign_end_id(ss, ss->ends()[1]);
    //std::cout << end_id2 << std::endl;

    return 1;
}

int
SecondaryStructureUnittest::test_to_str() {
    sstruct::SecondaryStructureFactory sf;
    auto ss = sf.get_structure("GGGG&CCCC", "((((&))))");
    //std::cout << ss->basepairs().size() << std::endl;
    auto s = ss->to_str();
    auto ss1 = sstruct::str_to_secondary_structure(s);
    if(ss->motifs("ALL").size() != ss1.motifs("ALL").size()) {
        return 0;
    }
    return 1;
}


int
SecondaryStructureUnittest::run() {
    if (test_creation_residue() == 0)  {  std::cout << "test_creation_residue failed" << std::endl; }
    if (test_creation() == 0)          {  std::cout << "test_creation failed" << std::endl; }
    if (test_get_residue() == 0)       {  std::cout << "test_get_residue failed" << std::endl; }
    if (test_assign_end_id() == 0)     {  std::cout << "test_assign_end_id failed" << std::endl; }
    if (test_to_str() == 0)            {  std::cout << "test_to_str failed" << std::endl; }
    return 1;
}
