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


//RNAMake Libraries
#include "base/types.h"

/*static
int
callback(void *NotUsed, int argc, char **argv, char **azColName){
    int i;
    for(i=0; i<argc; i++){
        printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
    }
    printf("\n");
    return 0;
}*/

class Sqlite3Connection {
public:
    Sqlite3Connection() {}
    
    Sqlite3Connection(String const &);
    
    ~Sqlite3Connection() {
        delete zErrMsg_;
        sqlite3_close(db_);
    }
    
public:
    
    void
    query(String const);
    
    int
    count();
    
    
public: //getters
    
    
protected:
    int rc_;
    int ic_;
    String query_statement_;
    std::string db_name_;
    sqlite3* db_;
    sqlite3_stmt* stmt_;
    char* zErrMsg_;
    int transaction_;

    
};

#endif /* defined(__RNAMake__sqlite3_connection__) */
