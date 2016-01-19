//
//  motif_sqlite_connection_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_sqlite_connection_unittest__
#define __RNAMake__motif_sqlite_connection_unittest__

#include <stdio.h>


#include "unittest.h"

class MotifSqliteConnectionUnittest : public Unittest {
public:
    MotifSqliteConnectionUnittest() {}
    
    ~MotifSqliteConnectionUnittest() {}
    
    
public:
    int
    test_creation();
    
    int
    test_next();
    
    void
    test_memory();
    
public:
    
    int
    run();
    
};


#endif /* defined(__RNAMake__motif_sqlite_connection_unittest__) */
