//
//  sqlite3_connection.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__segment_sqlite3_connection__
#define __RNAMake__segment_sqlite3_connection__

#include <stdio.h>
#include <iostream>
#include <sqlite3.h>

//RNAMake Libraries
#include "base/types.h"
#include "util/sqlite3_connection.h"

namespace resources {

struct SegmentSqliteData {
    SegmentSqliteData() :
            data(""), name(""), end_name(""), end_id(""),
            id("0") {}

    String data, name, end_name, end_id, id;

};

typedef std::shared_ptr<SegmentSqliteData> SegmentSqliteDataOP;

class SegmentSqliteConnection : public util::Sqlite3Connection {
public:
    SegmentSqliteConnection() {}

    SegmentSqliteConnection(
            String const & path):
            util::Sqlite3Connection(path),
            data_(std::make_shared<SegmentSqliteData>()) {}
public:

    SegmentSqliteDataOP const &
    next();

    SegmentSqliteDataOP const &
    contains();

    inline
    void
    clear() {
        if (rc_ == SQLITE_ROW || rc_ == SQLITE_DONE) {
            sqlite3_finalize(stmt_);
        }
    }

private:
    SegmentSqliteDataOP data_;
};

}

#endif /* defined(__RNAMake__segment_sqlite3_connection__) */
