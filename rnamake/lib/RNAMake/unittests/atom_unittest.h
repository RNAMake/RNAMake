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

class AtomUnittest : public Unittest {
public:
    AtomUnittest() {}
    
    ~AtomUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_to_pdb_str();
    
    int
    test_str_to_atom();
    
    int
    test_copy();
        
public:
    
    int
    run() {
        if (test_creation() == 0)    {  std::cout << "test_creation failed" << std::endl; }
        if (test_to_pdb_str() == 0)  {  std::cout << "test_to_pdb_str failed" << std::endl; }
        if (test_str_to_atom() == 0) {  std::cout << "test_str_to_atom failed" << std::endl; }
        if (test_copy() == 0)        {  std::cout << "test_copy failed" << std::endl; }
        
        return 1;
    }
    
    
};

#endif /* defined(__RNAMake__atom_unittest__) */
