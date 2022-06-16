//
//  sqlite3_connection.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_sqlite3_connection__
#define __RNAMake__motif_sqlite3_connection__

#include <iostream>
#include <sqlite3.h>
#include <stdio.h>

// RNAMake Libraries
#include "base/types.hpp"
#include "util/sqlite3_connection.h"

namespace resources {

struct MotifSqliteData {
  MotifSqliteData() : data(""), name(""), end_name(""), end_id(""), id("0") {}

  String data, name, end_name, end_id, id;
};

typedef std::shared_ptr<MotifSqliteData> MotifSqliteDataOP;

class MotifSqliteConnection : public util::Sqlite3Connection {
public:
  MotifSqliteConnection() {}

  MotifSqliteConnection(String const &path)
      : util::Sqlite3Connection(path),
        data_(std::make_shared<MotifSqliteData>()) {}

public:
  MotifSqliteDataOP const &next();

  MotifSqliteDataOP const &contains();

  inline void clear() {
    if (rc_ == SQLITE_ROW || rc_ == SQLITE_DONE) {
      sqlite3_finalize(stmt_);
    }
  }

private:
  MotifSqliteDataOP data_;
};

} // namespace resources

#endif /* defined(__RNAMake__motif_sqlite3_connection__) */
