//
// Created by Joseph Yesselman on 12/30/17.
//

#include <base/log.hpp>
#include <util/sqlite/connection.hpp>

namespace util::sqlite {

// interface to get rows //////////////////////////////////////////////////////
int Connection::setup_row_iteration(const String &command) const {
  if (!_done) {
    abort_iterate_rows();
  }
  _row = {};
  _first_row = true;
  _done = false;
  _prepare(command);
  if (sqlite3_step(_stmt) != SQLITE_ROW) {
    String msg = "cannot get sqlite row likely did not did not call setup "
                 "first ...";
    throw SqliteException(msg);
  }
  return 0;
}

const Row &Connection::next() const {
  const Row &r = _generate_row_from_statement();
  _rc = sqlite3_step(_stmt);
  if (_rc != SQLITE_ROW) {
    _done = true;
  }
  return r;
}

const Row &Connection::get_first_row(const String &command) const {
  _row = {};
  _first_row = true;
  _prepare(command);
  if (sqlite3_step(_stmt) != SQLITE_ROW) {
    String msg = "cannot get sqlite row with commmand  '" + command + "'";
    throw SqliteException(msg);
  }
  const Row &row = _generate_row_from_statement();
  sqlite3_finalize(_stmt);
  return row;
}

void Connection::abort_iterate_rows() const { sqlite3_finalize(_stmt); }

bool Connection::has_more_rows() const {
  if (_done) {
    return false;
  } else {
    return true;
  }
}

// send sql commands //////////////////////////////////////////////////////////
int Connection::exec(const String &command) const {
  // TODO Send error?
  _prepare(command);
  int err = _do_bind();
  if (err != SQLITE_OK) {
    sqlite3_finalize(_stmt);
    return err;
  }
  err = sqlite3_step(_stmt);
  if (err != SQLITE_DONE) {
    sqlite3_finalize(_stmt);
    return err;
  }
  return sqlite3_finalize(_stmt);
}

void Connection::start_transaction() const { _execute("BEGIN TRANSACTION;"); }

void Connection::commit_transaction() const { _execute("COMMIT;"); }

void Connection::rollback_transaction() const { _execute("ROLLBACK;"); }

// get table infomation ///////////////////////////////////////////////////////
TableDetailsOP Connection::get_table_details(const String &table_name) const {
  auto q = "SELECT sql FROM sqlite_master WHERE tbl_name = " +
           base::string::quoted(table_name);
  _prepare(q);
  _rc = sqlite3_step(_stmt);
  if (_rc != SQLITE_ROW) {
    throw SqliteException(table_name + " does not exist cannot get details");
  }
  auto table_str =
      String(reinterpret_cast<const char *>(sqlite3_column_text(_stmt, 0)));
  auto spl = base::string::split(table_str, "(");
  if (spl.size() != 3) {
    String msg = "not a valid table declaration: " + table_str;
    base::log_and_throw<SqliteException>(msg);
  }
  Strings columns = base::string::split(spl[1], ",");
  Strings names, types;
  std::map<String, bool> primary_keys;
  for (auto const &col : columns) {
    auto trimed_col = col;
    trimed_col = base::string::trim(trimed_col);
    auto col_spl = base::string::split(trimed_col, " ");
    if (col_spl.size() < 2) {
      throw SqliteException("cannot parse table declaration section: " + col);
    }

    if (col_spl[0] == "PRIMARY") {
      for (auto i = 1; i < col_spl.size(); i++) {
        primary_keys[col_spl[i]] = true;
      }
    } else {
      names.push_back(col_spl[0]);
      types.push_back(col_spl[1]);
    }
  }
  auto td = std::make_shared<TableDetails>(table_name);
  for (int i = 0; i < names.size(); i++) {
    auto is_primary = false;
    if (primary_keys.find(names[i]) == primary_keys.end()) {
      is_primary = true;
    }
    td->add_column(names[i], types[i], is_primary);
  }
  return td;
}

// binding ////////////////////////////////////////////////////////////////////
void Connection::bind(int param, const Blob &blob) {
  BindType b(param, &blob);
  _bind.push_back(b);
}

void Connection::bind(int param, const String &text) {
  BindType t(param, &text);
  _bind.push_back(t);
}

void Connection::bind(const char *param, const Blob &blob) {
  BindType b(param, &blob);
  _bind.push_back(b);
}

void Connection::bind(const char *param, const String &text) {
  BindType t(param, &text);
  _bind.push_back(t);
}

// private function ///////////////////////////////////////////////////////////
int Connection::_execute(const String &command) const {
  _rc = sqlite3_exec(_db(), command.c_str(), nullptr, nullptr, &_zErrMsg);
  if (_rc != SQLITE_OK) {
    throw SqliteException(String(_zErrMsg));
  }
  return _rc;
}

void Connection::_prepare(const String &command) const {
  _rc = sqlite3_prepare_v2(_db(), command.c_str(),
                           (int)strlen(command.c_str()) + 1, &_stmt, nullptr);
  if (_rc != SQLITE_OK) {
    throw SqliteException("error=" + std::to_string(_rc) +
                          " returned with query: " + command);
  }
}

const Row &Connection::_generate_row_from_statement() const {
  int col_type;
  if (_first_row) {
    _row.resize(sqlite3_column_count(_stmt));
  }
  int pos = 0;
  for (int i = 0; i < sqlite3_column_count(_stmt); ++i) {
    col_type = sqlite3_column_type(_stmt, i);
    if (col_type == SQLITE_INTEGER) {
      _row[pos] =
          Field(sqlite3_column_name(_stmt, i), sqlite3_column_int(_stmt, i));
    } else if (col_type == SQLITE_FLOAT) {
      _row[pos] =
          Field(sqlite3_column_name(_stmt, i), sqlite3_column_double(_stmt, i));
    } else if (col_type == SQLITE_BLOB) {
      auto const *blob = reinterpret_cast<const std::uint8_t *>(
          ::sqlite3_column_blob(_stmt, i));
      Blob b(&blob[0], &blob[::sqlite3_column_bytes(_stmt, i)]);
      _row[pos] = Field(sqlite3_column_name(_stmt, i), b);
    } else if (col_type == SQLITE3_TEXT) {
      _row[pos] =
          Field(sqlite3_column_name(_stmt, i),
                reinterpret_cast<const char *>(sqlite3_column_text(_stmt, i)));
    } else {
      throw SqliteException("not supported sqlite3_column type: " +
                            std::to_string(col_type));
    }
    pos++;
  }
  return _row;
}

int Connection::_do_bind() const {
  // TODO this should really be refactored uses a lot of old syntax
  int err = SQLITE_OK;
  for (auto &it : _bind) {
    if (it.type_ == SQLITE_BLOB) {
      const auto *v = static_cast<const std::vector<std::uint8_t> *>(it.ptr_);
      err = ::sqlite3_bind_blob(
          _stmt,
          !it.param_str_ ? it.param_num_
                         : ::sqlite3_bind_parameter_index(_stmt, it.param_str_),
          v->size() ? (const void *)&(*v)[0] : nullptr, (int)v->size(),
          SQLITE_STATIC);
      if (err != SQLITE_OK)
        break;
    }
    if (it.type_ == SQLITE_TEXT) {
      const auto *s = static_cast<const std::string *>(it.ptr_);
      err = ::sqlite3_bind_text(
          _stmt,
          !it.param_str_ ? it.param_num_
                         : ::sqlite3_bind_parameter_index(_stmt, it.param_str_),
          s->size() ? s->c_str() : nullptr, (int)s->size(), SQLITE_STATIC);
      if (err != SQLITE_OK)
        break;
    }
    if (it.type_ == SQLITE_NULL) {
      err = ::sqlite3_bind_null(
          _stmt, !it.param_str_
                     ? it.param_num_
                     : ::sqlite3_bind_parameter_index(_stmt, it.param_str_));
      if (err != SQLITE_OK)
        break;
    }
  }
  _bind.clear(); // clear bindings
  return err;
}

// wrapper functions ///////////////////////////////////////////////////////////
void create_table(const Connection &conn, const TableDetails &td) {
  auto table_str = "CREATE TABLE " + td.name() + "(";
  auto i = 0;
  for (auto const &col : td) {
    table_str += col.name + " " + col.type;
    i++;
    if (i != td.size()) {
      table_str += ", ";
    }
  }
  if (td.has_primary_key()) {
    table_str += ", PRIMARY KEY( ";
    for (auto const &col : td) {
      if (col.is_primary) {
        table_str += col.name + " ";
      }
    }
    table_str += ") ";
  }
  table_str += ")";
  conn.exec(table_str);
  LOG_INFO << "creating table: " + table_str +
                  " in database: " + conn.get_database_name();
}

void insert_many(const Connection &conn, const String &table_name,
                 const std::vector<Strings> &data) {

  TableDetailsOP td = conn.get_table_details(table_name);
  String insert_str = "INSERT INTO " + table_name + "( ";
  int i = 0;
  for (auto const &col : *td) {
    insert_str += col.name;
    i++;
    if (i != td->size()) {
      insert_str += ", ";
    }
  }
  insert_str += ") VALUES ";
  conn.start_transaction();
  String insert_statement;
  for (auto const &row : data) {
    i = -1;
    insert_statement = insert_str + "(";
    for (auto const &element : row) {
      i++;
      insert_statement += "'" + element + "'";
      if (i != row.size() - 1) {
        insert_statement += ",";
      }
    }
    insert_statement += ");";
    conn.exec(insert_statement);
  }
  conn.commit_transaction();
}

} // namespace util::sqlite