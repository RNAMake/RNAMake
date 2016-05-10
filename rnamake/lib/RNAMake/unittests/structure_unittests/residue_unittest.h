//
//  residue_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__residue_unittest__
#define __RNAMake__residue_unittest__

#include <stdio.h>
#include <iostream>

//RNAMake Headers
#include "structure/residue.h"
#include "structure/residue_type_set.h"

#include "unittest.h"

namespace unittests {

class ResidueUnittest : public Unittest {
public:
    ResidueUnittest();

    ~ResidueUnittest() {}
    
public:
    
    int
    size() { return 7; }
    
    void
    test_bead_creation();
    
    void
    test_get_atom();
    
    void
    test_connected_to();
    
    void
    test_get_beads();
    
    void
    test_copy();
    
    void
    test_to_str();
    
    void
    test_equals();
    
    void
    test_memory_management();
    
public:
    
    int
    run();
    
    int
    run_all();
    
private:
    
    ResidueTypeSet rts_;
    ResidueOPs residues_;
    
    
};

    
}

#endif /* defined(__RNAMake__residue_unittest__) */
