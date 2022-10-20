//
//  sqlite3_connection.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__segment_sqlite3_connection__
#define __RNAMake__segment_sqlite3_connection__

#include <iostream>
#include <sqlite3.h>
#include <stdio.h>

// RNAMake Libraries
#include "base/types.hpp"
#include "util/sqlite/connection.hpp"

namespace resources {

struct SegmentSqliteData {
  SegmentSqliteData() : data(""), name(""), end_name(""), end_id(""), id("0") {}

  String data, name, end_name, end_id, id;
};

typedef std::shared_ptr<SegmentSqliteData> SegmentSqliteDataOP;

class SegmentSqliteConnection : public util::Sqlite3Connection {
public:
  SegmentSqliteConnection() {}

  SegmentSqliteConnection(String const &path)
      : util::Sqlite3Connection(path),
        _data(std::make_shared<SegmentSqliteData>()) {}

public:
  SegmentSqliteDataOP const &next();

  SegmentSqliteDataOP const &contains();

  inline void clear() {
    if (_rc == SQLITE_ROW || _rc == SQLITE_DONE) {
      sqlite3_finalize(_stmt);
    }
  }

private:
  SegmentSqliteDataOP _data;
};

} // namespace resources

#endif /* defined(__RNAMake__segment_sqlite3_connection__) */
