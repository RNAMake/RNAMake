//
//  sqlite3_connection.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sqlite3_connection__
#define __RNAMake__sqlite3_connection__

#include <stdio.h>
#include <iostream>
#include <sqlite3.h>

namespace util {

//RNAMake Libraries
#include "base/types.h"

class Sqlite3ConnectionException : public std::runtime_error {
public:
    Sqlite3ConnectionException(String const & message) :
            std::runtime_error(message) {}
};


class Sqlite3Connection {
public:
    Sqlite3Connection() :
            setup_(0) {}

    Sqlite3Connection(String const &);

    ~Sqlite3Connection() {
        if (setup_) {
            delete zErrMsg_;
            sqlite3_close(db_);
        }
    }

public:

    void
    query(String const &);

    int
    count();

    Strings
    fetch_one(String const &);


public: //getters


protected:
    int rc_;
    int ic_;
    int setup_;
    String query_statement_;
    String db_name_;
    sqlite3 *db_;
    sqlite3_stmt *stmt_;
    char *zErrMsg_;
};

}

#endif /* defined(__RNAMake__sqlite3_connection__) */
