//
//  secondary_structure_factory_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure_factory_unittest.h"
#include "secondary_structure/secondary_structure_factory.h"

namespace unittests {

void
SecondaryStructureFactoryUnittest::test_creation() {
    sstruct::SecondaryStructureFactory ssf;
    auto m = ssf.motif("GGGGG+CCCCC", "(((((+)))))");
    if(m->residues().size() != 10) {
        throw UnittestException("did not get the expected number of res");
    }
    
    if(m->get_residue(10, "B", "") == nullptr) {
        throw UnittestException("could not find last residue");
    }
    
    
    /*if(ss->motifs("BP_STEP").size() != 4 || ss->motifs("ALL").size() != 4) {
        return 0;
    }

    auto ss1 = sf.get_structure("GGGAAGG+CCCCC", "(((..((+)))))");
    if(ss1->motifs("BULGE").size() != 1) {
        return 0;
    }
    
    auto ss2 = ss1->copy();
    if(ss2.motifs("BP_STEP").size() != 3 || ss2.motifs("ALL").size() != 4 ||
       ss2.motifs("BULGE").size() != 1 ) {
        return 0;
    }*/
}


int
SecondaryStructureFactoryUnittest::run() {
    test_creation();
    return 1;
    
}

}