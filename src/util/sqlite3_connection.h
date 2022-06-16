//
//  sqlite3_connection.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sqlite3_connection__
#define __RNAMake__sqlite3_connection__

#include <iostream>
#include <sqlite3.h>
#include <stdexcept>
#include <stdio.h>
// RNAMake Libraries
#include "base/types.hpp"

namespace util {

class Sqlite3ConnectionException : public std::runtime_error {
public:
  Sqlite3ConnectionException(String const &message)
      : std::runtime_error(message) {}
};

class Sqlite3Connection {
public:
  Sqlite3Connection() : _setup(0) {}

  Sqlite3Connection(String const &);

  ~Sqlite3Connection() {
    if (_setup) {
      delete _zErrMsg;
      sqlite3_close(_db);
    }
  }

public:
  void query(String const &);

  int count();

  Strings fetch_one(String const &);

public: // getters
protected:
  int _rc;
  int _ic;
  int _setup;
  String _query_statement;
  String _db_name;
  sqlite3 *_db;
  sqlite3_stmt *_stmt;
  char *_zErrMsg;
};

} // namespace util

#endif /* defined(__RNAMake__sqlite3_connection__) */
