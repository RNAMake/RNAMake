//
//  x3dna_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__x3dna_unittest__
#define __RNAMake__x3dna_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"
#include "util/x3dna.h"

class X3dnaUnittest : public Unittest {
public:
    X3dnaUnittest() {  x_ = X3dna();
    }
    
    ~X3dnaUnittest() {} 

public:
    
    int
    test_creation();

    int
    test_get_ref_frame();
    
    int
    test_generate_ref_frame();
    
    int
    test_generate_dssr_file();
    
    int
    test_res_compare();
    
    int
    test_get_basepairs();
    
public:
    
    int
    run();
    
private:
    X3dna x_;
    
};


#endif /* defined(__RNAMake__x3dna_unittest__) */
