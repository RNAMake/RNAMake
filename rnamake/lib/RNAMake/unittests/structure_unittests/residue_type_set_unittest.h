//
//  residue_type_set_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__residue_type_set_unittest__
#define __RNAMake__residue_type_set_unittest__

#include <stdio.h>
#include <iostream>

//RNAMake Headers
#include "base/types.h"
#include "structure/residue_type_set.h"
#include "structure/residue_type.h"

#include "unittest.h"


class ResidueTypeSetUnittest : public Unittest {
public:
    ResidueTypeSetUnittest() {}
    
    ~ResidueTypeSetUnittest() {}
    
public:
    
    int
    test_creation_residue_type();
    
    int
    test_match_name();
    
    int
    test_creation();
    
    int
    test_get_rtype_by_resname();
    
    int
    test_atom_pos_by_name();
    
public:
    
    int
    run() {
        if (test_creation_residue_type() == 0) { std::cout << "test_creation_residue_type failed" << std::endl;  }
        if (test_creation() == 0 )             { std::cout << "test_creation failed" << std::endl; }
        if (test_match_name() == 0 )           { std::cout << "test_match_name failed" << std::endl; }
        if (test_get_rtype_by_resname() == 0 ) { std::cout << "test_get_rtype_by_resname failed" << std::endl; }
        if (test_atom_pos_by_name() == 0 )     { std::cout << "test_atom_pos_by_name failed" << std::endl; }
        
        return 1;
    }
    
    
};


#endif /* defined(__RNAMake__residue_type_set_unittest__) */
