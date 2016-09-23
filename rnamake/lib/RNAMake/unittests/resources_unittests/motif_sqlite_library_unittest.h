//
//  motif_sqlite_library_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_sqlite_library_unittest__
#define __RNAMake__motif_sqlite_library_unittest__

#include <stdio.h>

#include "unittest.h"
#include "resources/motif_sqlite_library.h"

class MotifSqliteLibraryUnittest : public Unittest {
public:
    
    MotifSqliteLibraryUnittest() {}
    
    ~MotifSqliteLibraryUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_get();
    
    int
    test_get_random();
    
    int
    test_all();
    
    int
    test_get_multi();
    
    int
    test_contains();
    
    void
    test_memory();
    
public:
    
    int
    run();
    
    
};

#endif /* defined(__RNAMake__motif_sqlite_library_unittest__) */
