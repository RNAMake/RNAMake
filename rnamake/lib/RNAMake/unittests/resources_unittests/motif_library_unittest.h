//
//  motif_library_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_library_unittest__
#define __RNAMake__motif_library_unittest__

#include <stdio.h>

#include "unittest.h"
#include "resources/motif_library.h"

class MotifLibraryUnittest : public Unittest {
public:
    
    MotifLibraryUnittest() {}
    
    ~MotifLibraryUnittest() {}
    
public:
    
    int
    test_load_all();
    
    int
    test_get_motif();
    
    int
    test_contains_motif();
    
public:
    
    int
    run();

    
};

#endif /* defined(__RNAMake__motif_library_unittest__) */
