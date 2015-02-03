//
//  sqlite3_connection.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "sqlite3_connection.h"
#include "settings.h"

Sqlite3Connection::Sqlite3Connection(
    MotifType const & mtype) {
    
    zErrMsg_ = 0;
    ic_ = 0;
    String path = get_db_path(mtype);
    rc_ = sqlite3_open(path.c_str(), &db_);
    if( rc_ ){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_));
        exit(0);
    }
    values_ = Strings(2);

}

String
Sqlite3Connection::get_db_path(MotifType const & mtype) {
    String path = resources_path() + "/motif_libraries/";
    if     ( mtype == TWOWAY) { path += "two_ways"; }
    else if( mtype == HELIX ) { path += "helices";  }
    else { throw "could not find mtype"; }
    return path+".db";
    
}

int
Sqlite3Connection::next() {
    if(rc_ != SQLITE_ROW) {
        sqlite3_finalize(stmt_);
        return 0;
    }
    values_[0] = String(reinterpret_cast<const char*>(sqlite3_column_text(stmt_,0)));
    values_[1] = String(reinterpret_cast<const char*>(sqlite3_column_text(stmt_,1)));

    return 1;
}

void
Sqlite3Connection::query(String const query_statement) {
    rc_ = sqlite3_prepare_v2(db_,
                             query_statement.c_str(),
                             strlen(query_statement.c_str())+1,
                             &stmt_,
                             NULL);
    rc_ = sqlite3_step(stmt_);
}


