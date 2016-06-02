//
//  sqlite3_connection.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <string.h>
#include "base/settings.h"
#include "util/sqlite3_connection.h"


Sqlite3Connection::Sqlite3Connection(
    String const & path) {
    
    zErrMsg_ = 0;
    ic_ = 0;
    rc_ = sqlite3_open(path.c_str(), &db_);
    if( rc_ ){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_));
        exit(0);
    }

}

void
Sqlite3Connection::query(String const & query_statement) {
    
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


