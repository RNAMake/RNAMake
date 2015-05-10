//
//  sqlite3_connection.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_sqlite3_connection__
#define __RNAMake__motif_sqlite3_connection__

#include <stdio.h>
#include <iostream>
#include <sqlite3.h>

//RNAMake Libraries
#include "base/types.h"
#include "util/sqlite3_connection.h"

class MotifSqliteConnection : public Sqlite3Connection {
public:
    MotifSqliteConnection(String const & path):
    Sqlite3Connection(path)
    {}
    
    ~MotifSqliteConnection() {
        delete zErrMsg_;
        sqlite3_close(db_);
    }
    
public:
    
    Strings const &
    next();
    
};

#endif /* defined(__RNAMake__motif_sqlite3_connection__) */
