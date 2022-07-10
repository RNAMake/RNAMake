//
// Created by Joseph Yesselman on 1/3/18.
//

#include <resource_management/sqlite_library.h>

namespace resource_management {

SqliteLibrary::SqliteLibrary(const String &db_path, const String &table_name):
      _conn(util::sqlite::Connection(db_path)),
      _table_details(util::sqlite::TableDetails(table_name)) {

  if (_conn.is_database_created()) {
    std::remove(db_path.c_str());
    throw ResourceManagementException(
        "SqliteLibrary expects an already built sqlite3 data with a given "
        "table, was given db_path: " +
        db_path + " which does not exist!");
  }

  _table_details = *(_conn.get_table_details(table_name));
}

size_t SqliteLibrary::get_num_of_rows() {
  _query_string = "SELECT COUNT(*) FROM " + _table_details.name();
  auto &row = _conn.get_first_row(_query_string);
  int count = row[0].get_int();
  return (size_t)count;
}

void SqliteLibrary::_generate_query(
    const Strings &retrieved_columns,
    const StringStringMap &restraint_col_and_vals) const {

  _query_string = "SELECT ";
  int i = 0;
  for (auto const &col_name : retrieved_columns) {
    i++;
    if (!_is_valid_key(col_name)) {
      throw ResourceManagementException(
          "attempting to use colunn key: " + col_name +
          " does not exist in sqlite library");
    }
    _query_string += col_name + " ";
    if (i != retrieved_columns.size()) {
      _query_string += ",";
    }
  }
  _query_string += " FROM " + _table_details.name() + " WHERE ";
  i = 0;
  for (auto const &kv : restraint_col_and_vals) {
    i++;
    if (!_is_valid_key(kv.first)) {
      throw ResourceManagementException(
          "attempting to use colunn key: " + kv.first +
          " does not exist in sqlite library");
    }
    _query_string += kv.first + "='" + kv.second + "' ";
    if (i != restraint_col_and_vals.size()) {
      _query_string += "AND ";
    }
  }
}

bool SqliteLibrary::_is_valid_key(const String &key) const {
  for (auto const &col : _table_details) {
    if (col.name == key) {
      return true;
    }
  }
  return false;
}

bool SqliteLibrary::_does_query_return_rows(
    const StringStringMap &restraint_col_and_vals) const {
  _query_string = "SELECT COUNT(*) FROM " + _table_details.name() + " WHERE ";
  int i = 0;
  for (auto const &kv : restraint_col_and_vals) {
    i++;
    //if (!_is_valid_key(kv.first)) {
    //  throw ResourceManagementException(
    //      "attempting to use colunn key: " + kv.first +
    //      " does not exist in sqlite library");
    //}
    _query_string += kv.first + "='" + kv.second + "' ";
    if (i != restraint_col_and_vals.size()) {
      _query_string += "AND ";
    }
  }
  const auto &row = _conn.get_first_row(_query_string);
  int count = row[0].get_int();
  if (count == 0) {
    return false;
  } else {
    return true;
  }
}

} // namespace resource_management