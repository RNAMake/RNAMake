//
// Created by Joseph Yesselman on 12/30/17.
//

#ifndef RNAMAKE_NEW_SQLITE_DATABASE_H
#define RNAMAKE_NEW_SQLITE_DATABASE_H

#include <sqlite3.h>

#include <base/paths.hpp>
#include <base/types.hpp>

namespace util::sqlite {

class Database {
public:// construction
  inline explicit Database(const String &name) { _connect(name); }

  inline ~Database() {
    // close the db
    (void) close();
  }

  // cannot copy
  Database &operator=(const Database &) = delete;


public:
  // close the database
  int close();

public:// getters
  /// @brief Bool is the database file open?
  [[nodiscard]] inline bool is_open() const { return _open; }

  /// @brief Bool did the database need to be created
  [[nodiscard]] inline bool is_created() const { return _created; }

  // SQLite3 access
  [[nodiscard]] inline const sqlite3 *get() const { return _db; }

  // @brief get name of database file
  [[nodiscard]] inline const String &get_name() const { return _name; }

public:
  inline sqlite3 *operator()() const { return _db; }

private:
  // open (connect) the database
  void _connect(const String &name,
                int flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE);


private:
  sqlite3 *_db = nullptr;// associated db
  String _name;          // db filename
  bool _open = false;    // db open status
  bool _created = false;
};

}// namespace util::sqlite

#endif//RNAMAKE_NEW_SQLITE_DATABASE_H
