//
//  sqlite3_connection.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//


#include "util/settings.h"
#include "resources/motif_sqlite_connection.h"


Strings const &
MotifSqliteConnection::next() {
    if(rc_ != SQLITE_ROW) {
        sqlite3_finalize(stmt_);
        values_[0] = "";
        values_[1] = "";
        return values_;
    }
    values_[0] = String(reinterpret_cast<const char*>(sqlite3_column_text(stmt_,0)));
    values_[1] = String(reinterpret_cast<const char*>(sqlite3_column_text(stmt_,1)));
    rc_ = sqlite3_step(stmt_);
    
    return values_;
}