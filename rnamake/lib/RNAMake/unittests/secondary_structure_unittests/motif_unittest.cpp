//
//  motif_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_unittest.h"
#include "secondary_structure/secondary_structure_factory.h"

namespace unittests {
namespace sstruct_unittests  {

void
MotifUnittest::test_creation() {
    sstruct::SecondaryStructureFactory ssf;
    auto m = ssf.motif("GGGGG+CCCCC", "(((((+)))))");
}
    
void
MotifUnittest::test_copy() {
    sstruct::SecondaryStructureFactory ssf;
    auto m = ssf.motif("GGGGG+CCCCC", "(((((+)))))");
    auto m_copy = sstruct::Motif(*m);
    for(auto const & r : m->residues()) {
        if(m_copy.get_residue(r->uuid()) == nullptr) {
            throw UnittestException("cannot find residue in copy");
        }
    }
    
    if(m->basepairs().size() != m_copy.basepairs().size()) {
        throw UnittestException("did not get the correct number of baseapirs in copy");
    }
    
    if(m->ends().size() != m_copy.ends().size()) {
        throw UnittestException("did not get the correct number of ends in copy");
    }
    
    //auto m1 = std::make_unique<const sstruct::Motif>(m_copy);
    //std::unique_ptr<const sstruct::Motif> & m2 = m1;
    //m2->mtype(HELIX);
    //std::cout << type_to_str(m1->mtype()) << std::endl;
}
    
void
MotifUnittest::test_to_str() {
    sstruct::SecondaryStructureFactory ssf;
    auto m = ssf.motif("GGGGG+CCCCC", "(((((+)))))");
    auto s = m->to_str();
    auto m1 = sstruct::Motif(*s);
    for(auto const & r : m->residues()) {
        if(m1.get_residue(r->num(), r->chain_id(), r->i_code()) == nullptr) {
            throw UnittestException("cannot find residue in copy");
        }
    }
    
}
    
    

int
MotifUnittest::run() {
    test_creation();
    test_copy();
    test_to_str();
    return 0;
}
    
} // sstruct_unittests
} // unittests























