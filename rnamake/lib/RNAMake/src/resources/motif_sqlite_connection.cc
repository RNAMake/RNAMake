//
//  sqlite3_connection.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//


#include "util/settings.h"
#include "resources/motif_sqlite_connection.h"


MotifSqliteDataOP const &
MotifSqliteConnection::next() {
    if(rc_ != SQLITE_ROW) {
        sqlite3_finalize(stmt_);
        data_->data = "";
        return data_;
    }
    data_->data     = String(reinterpret_cast<const char*>(sqlite3_column_text(stmt_,0)));
    data_->name     = String(reinterpret_cast<const char*>(sqlite3_column_text(stmt_,1)));
    data_->end_name = String(reinterpret_cast<const char*>(sqlite3_column_text(stmt_,2)));
    data_->end_id   = String(reinterpret_cast<const char*>(sqlite3_column_text(stmt_,3)));
    data_->id       = String(reinterpret_cast<const char*>(sqlite3_column_text(stmt_,4)));
    //values_[0] = String(reinterpret_cast<const char*>(sqlite3_column_text(stmt_,0)));
    //values_[1] = String(reinterpret_cast<const char*>(sqlite3_column_text(stmt_,1)));
    rc_ = sqlite3_step(stmt_);
    
    return data_;
}