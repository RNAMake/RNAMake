//
//  library_manager_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__library_manager_unittest__
#define __RNAMake__library_manager_unittest__

#include <stdio.h>

#include "unittest.h"

class LibraryManagerUnittest : public Unittest {
public:
    LibraryManagerUnittest();

    ~LibraryManagerUnittest() {}

    
public:
    int
    test_get_motif();
    
    int
    test_add_motif();
    
public:
    
    int
    run();
    
    void
    run_all();
    
    
private:
    
};
    
#endif /* defined(__RNAMake__library_manager_unittest__) */