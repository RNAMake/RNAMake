//
//  sqlite3_connection_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sqlite3_connection_unittest__
#define __RNAMake__sqlite3_connection_unittest__

#include <stdio.h>


//RNAMake Headers
#include "unittest.h"

class Sqlite3ConnectionUnittest : public Unittest {
public:
    Sqlite3ConnectionUnittest() {}
    
    ~Sqlite3ConnectionUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_count();
    
public:
    
    int
    run();
    
};

#endif /* defined(__RNAMake__sqlite3_connection_unittest__) */
