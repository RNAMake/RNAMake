//
//  added_motif_library_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/9/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__added_motif_library_unittest__
#define __RNAMake__added_motif_library_unittest__

#include <stdio.h>

#include "unittest.h"


class AddedMotifLibraryUnittest : public Unittest {
public:
    
    AddedMotifLibraryUnittest() {}
    
    ~AddedMotifLibraryUnittest() {}
    
public:
    
    int
    test_add_motif();
    
    int
    test_get();
    
    int
    test_get_multi();
    
    int
    test_contains();
    
public:
    
    int
    run();
    
};


#endif /* defined(__RNAMake__added_motif_library_unittest__) */
