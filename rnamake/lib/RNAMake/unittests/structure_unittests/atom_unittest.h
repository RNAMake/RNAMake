//
//  atom_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__atom_unittest__
#define __RNAMake__atom_unittest__

#include <stdio.h>
#include <iostream>

//RNAMake Headers
#include "unittest.h"
#include "math/xyz_vector.h"
#include "math/numerical.h"
#include "structure/atom.h"

namespace unittests {

class AtomUnittest : public Unittest {
public:

    AtomUnittest() {}
    
    ~AtomUnittest() {}
    
public:
    
    int
    size() { return 5; }
    
    int
    test_creation();
    
    int
    test_to_str();
    
    int
    test_to_pdb_str();
    
    int
    test_str_to_atom();
    
    int
    test_copy();
    
        
public:
    
    int
    run();
    
    int
    run_all();
    
    
};

}

#endif /* defined(__RNAMake__atom_unittest__) */
