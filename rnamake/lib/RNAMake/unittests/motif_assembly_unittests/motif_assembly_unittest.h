//
//  motif_assembly_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_assembly_unittest__
#define __RNAMake__motif_assembly_unittest__

#include <stdio.h>

#include "unittest.h"

class MotifAssemblyUnittest : public Unittest {
public:
    MotifAssemblyUnittest() {}
    
    ~MotifAssemblyUnittest() {}
    
public:
    
    int
    test_add_motif();
    
    int
    test_motif_end_indentity();
    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};

#endif /* defined(__RNAMake__motif_assembly_unittest__) */
