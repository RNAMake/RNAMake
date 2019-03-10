//
//  sqlite3_connection.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <string.h>
#include "base/settings.h"
#include "base/file_io.h"
#include "util/sqlite3_connection.h"


Sqlite3Connection::Sqlite3Connection(
    String const & path) {
    
    zErrMsg_ = 0;
    ic_ = 0;
    if(!base::file_exists(path)) {
        throw Sqlite3ConnectionException(
            "Can't open sqlite3 database " + path + " it does not exist");
    }
    
    rc_ = sqlite3_open(path.c_str(), &db_);
    if( rc_ ){
        throw Sqlite3ConnectionException("Can't open sqlite3 database " + path);
    }
    
    setup_ = 1;

}

void
Sqlite3Connection::query(String const & query_statement) {
    if(setup_ == 0) {
        throw Sqlite3ConnectionException("Cannot call query if no database file is supplied!");
    }
    
    rc_ = sqlite3_prepare_v2(db_,
                             query_statement.c_str(),
                             (int)strlen(query_statement.c_str())+1,
                             &stmt_,
                             NULL);
    rc_ = sqlite3_step(stmt_);

}

int
Sqlite3Connection::count() {
    query("SELECT count(*) from data_table");
    int count = sqlite3_column_int(stmt_,0);
    rc_ = sqlite3_step(stmt_);
    return count;
}



Strings
Sqlite3Connection::fetch_one(String const & query_statement) {
    query(query_statement);
    
    auto results = Strings();
    
    int col = sqlite3_column_count(stmt_);
    
    for(int i = 0; i < col; i++) {
        results.push_back(String(reinterpret_cast<const char*>(sqlite3_column_text(stmt_,i))));
    }
    
    sqlite3_finalize(stmt_);
    return results;
    
    
}

