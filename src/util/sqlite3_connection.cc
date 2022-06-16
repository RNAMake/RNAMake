//
//  sqlite3_connection.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <string.h>
//#include "base/settings.h"
//#include "base/file_io.h"
#include "util/sqlite3_connection.h"

namespace util {

Sqlite3Connection::Sqlite3Connection(String const &path) {

  _zErrMsg = 0;
  _ic = 0;

  if (!base::file_exists(path)) {
    throw Sqlite3ConnectionException("Can't open sqlite3 database " + path +
                                     " it does not exist");
  }

  _rc = sqlite3_open(path.c_str(), &_db);
  if (_rc) {
    throw Sqlite3ConnectionException("Can't open sqlite3 database " + path);
  }

  _setup = 1;
}

void Sqlite3Connection::query(String const &query_statement) {
  if (_setup == 0) {
    throw Sqlite3ConnectionException(
        "Cannot call query if no database file is supplied!");
  }

  _rc = sqlite3_prepare_v2(_db, query_statement.c_str(),
                           (int)strlen(query_statement.c_str()) + 1, &_stmt,
                           NULL);
  _rc = sqlite3_step(_stmt);
}

int Sqlite3Connection::count() {
  query("SELECT count(*) from data_table");
  int count = sqlite3_column_int(_stmt, 0);
  _rc = sqlite3_step(_stmt);
  return count;
}

Strings Sqlite3Connection::fetch_one(String const &query_statement) {
  query(query_statement);

  auto results = Strings();

  int col = sqlite3_column_count(_stmt);

  for (int i = 0; i < col; i++) {
    results.push_back(
        String(reinterpret_cast<const char *>(sqlite3_column_text(_stmt, i))));
  }

  sqlite3_finalize(_stmt);
  return results;
}
} // namespace util
