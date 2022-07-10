//
// Created by Joseph Yesselman on 1/3/18.
//

#ifndef RNAMAKE_NEW_SQLITE_LIBRARY_H
#define RNAMAKE_NEW_SQLITE_LIBRARY_H

#include <resource_management/exception.hpp>
#include <util/sqlite/connection.hpp>

namespace resource_management {

class SqliteLibrary {
public:
  SqliteLibrary(const String &, const String &);

  virtual ~SqliteLibrary() = default;

public:
  size_t get_num_of_rows();

protected:
  void _generate_query(const Strings &, const StringStringMap &) const;

  bool _is_valid_key(const String &) const;

  bool _does_query_return_rows(const StringStringMap &) const;

protected:
  util::sqlite::TableDetails _table_details;
  // sqlite3 api commands cannot be const
  mutable util::sqlite::Connection _conn;
  mutable String _query_string;
};

} // namespace resource_management

#endif // RNAMAKE_NEW_SQLITE_LIBRARY_H
