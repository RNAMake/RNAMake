//
// Created by Joseph Yesselman on 12/30/17.
//

#ifndef RNAMAKE_NEW_SQLITE3_CONNECTION_H
#define RNAMAKE_NEW_SQLITE3_CONNECTION_H

#include <util/sqlite/database.h>
#include <util/sqlite/field.h>
#include <util/sqlite/table_details.h>

namespace util::sqlite {

class Connection {
public:
  inline explicit Connection(const Database &db) : db_(db) {}

  ~Connection() = default;

public:// interface to get rows
  int setup_row_iteration(const String &);

  RowOP next() {
    // TODO throw error here!!
    if (sqlite3_step(stmt_) != SQLITE_ROW) { return RowOP(nullptr); }
    return _generate_row_from_statement();
  }


  RowOP get_first_row(const String &command) {
    int err = _prepare(command);
    if (sqlite3_step(stmt_) != SQLITE_ROW) {
      sqlite3_finalize(stmt_);
      return {nullptr};
    }
    auto row = _generate_row_from_statement();
    sqlite3_finalize(stmt_);
    return row;
  }

  void abort_iterate_rows() { sqlite3_finalize(stmt_); }

public:// get other table infos
  TableDetailsOP get_table_details(String const &table_name) {

    auto q = "SELECT sql FROM sqlite_master WHERE tbl_name = " +
             _quoted_string(table_name);
    _prepare(q);
    rc_ = sqlite3_step(stmt_);
    if (rc_ != SQLITE_ROW) {
      throw SqliteException(table_name + " does not exist cannot get details");
    }

    auto table_str = String(
            reinterpret_cast<const char *>(sqlite3_column_text(stmt_, 0)));
    auto spl = base::string::split(table_str, "(");
    if (spl.size() != 3) {
      throw SqliteException("not a valid table declaration: " + table_str);
    }
    auto columns = base::string::split(spl[1], ",");
    auto names = Strings();
    auto types = Strings();
    auto primary_keys = std::map<String, bool>();
    for (auto const &col: columns) {
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

  inline String const &get_database_name() const { return db_.get_name(); }

  inline bool is_database_created() const { return db_.is_created(); }


public:
  int exec(String const &command) {
    int err = _prepare(command);
    err = do_bind();
    if (err != SQLITE_OK) {
      sqlite3_finalize(stmt_);
      return err;
    }
    err = sqlite3_step(stmt_);
    if (err != SQLITE_DONE) {
      sqlite3_finalize(stmt_);
      return err;
    }
    return sqlite3_finalize(stmt_);
  }

  void start_transaction() { _execute("BEGIN TRANSACTION;"); }

  void commit_transaction() { _execute("COMMIT;"); }

  void rollback_transaction() { _execute("ROLLBACK;"); }

  // bind a BLOB or TEXT to query
  // CAUTION: vector and string MUST BE constant until end of query execution of exec()/use()/store()
  void bind(int param, std::vector<std::uint8_t> const &blob) {
    BindType b(param, &blob);
    bind_.push_back(b);
  }
  void bind(int param, String const &text) {
    BindType t(param, &text);
    bind_.push_back(t);
  }
  void bind(const char *param, std::vector<std::uint8_t> const &blob) {
    BindType b(param, &blob);
    bind_.push_back(b);
  }
  void bind(const char *param, std::string const &text) {
    BindType t(param, &text);
    bind_.push_back(t);
  }

private:
  int _execute(String const &command) {
    rc_ = sqlite3_exec(db_(), command.c_str(), nullptr, 0, &zErrMsg_);
    if (rc_ != SQLITE_OK) { throw SqliteException(String(zErrMsg_)); }
    return rc_;
  }

  int _prepare(String const &command) {
    rc_ = sqlite3_prepare_v2(db_(), command.c_str(),
                             (int) strlen(command.c_str()) + 1, &stmt_, NULL);
    if (rc_ != SQLITE_OK) {
      throw SqliteException("error=" + std::to_string(rc_) +
                            " returned with query: " + command);
    }
    return rc_;
  }

  RowOP _generate_row_from_statement() {
    int col_type;
    if (_first_row) { _row.resize(sqlite3_column_count(stmt_)); }
    int pos = 0;
    for (int i = 0; i < sqlite3_column_count(stmt_); ++i) {
      col_type = sqlite3_column_type(stmt_, i);
      if (col_type == SQLITE_INTEGER) {
        _row[pos] = Field(sqlite3_column_name(stmt_, i),
                          sqlite3_column_int(stmt_, i));
      } else if (col_type == SQLITE_FLOAT) {
        _row[pos] = Field(sqlite3_column_name(stmt_, i),
                          sqlite3_column_double(stmt_, i));
      } else if (col_type == SQLITE_BLOB) {
        auto const *blob = reinterpret_cast<const std::uint8_t *>(
                ::sqlite3_column_blob(stmt_, i));
        std::vector<std::uint8_t> v(&blob[0],
                                    &blob[::sqlite3_column_bytes(stmt_, i)]);
        _row[pos] = Field(sqlite3_column_name(stmt_, i), v);
      } else if (col_type == SQLITE3_TEXT) {
        _row[pos] = Field(
                sqlite3_column_name(stmt_, i),
                reinterpret_cast<const char *>(sqlite3_column_text(stmt_, i)));
      } else {
        throw SqliteException("not supported sqlite3_column type: " +
                              std::to_string(col_type));
      }
      pos++;
    }
    return std::make_shared<Row>(_row);
  }


  String _quoted_string(String const &s) {
    return String("'") + s + String("'");
  }

  int do_bind() {
    int err = SQLITE_OK;
    for (std::vector<BindType>::const_iterator it = bind_.begin();
         it != bind_.end(); ++it) {
      if (it->type_ == SQLITE_BLOB) {
        const std::vector<std::uint8_t> *v =
                static_cast<const std::vector<std::uint8_t> *>(it->ptr_);
        err = ::sqlite3_bind_blob(
                stmt_,
                !it->param_str_
                        ? it->param_num_
                        : ::sqlite3_bind_parameter_index(stmt_, it->param_str_),
                v->size() ? (const void *) &(*v)[0] : nullptr, (int) v->size(),
                SQLITE_STATIC);
        if (err != SQLITE_OK) break;
      }
      if (it->type_ == SQLITE_TEXT) {
        const std::string *s = static_cast<const std::string *>(it->ptr_);
        err = ::sqlite3_bind_text(
                stmt_,
                !it->param_str_
                        ? it->param_num_
                        : ::sqlite3_bind_parameter_index(stmt_, it->param_str_),
                s->size() ? s->c_str() : nullptr, (int) s->size(),
                SQLITE_STATIC);
        if (err != SQLITE_OK) break;
      }
      if (it->type_ == SQLITE_NULL) {
        err = ::sqlite3_bind_null(
                stmt_, !it->param_str_ ? it->param_num_
                                       : ::sqlite3_bind_parameter_index(
                                                 stmt_, it->param_str_));
        if (err != SQLITE_OK) break;
      }
    }
    bind_.clear();// clear bindings
    return err;
  }


private:
  struct BindType {
    void const *ptr_;
    char const *param_str_;
    int param_num_;
    int type_;

    BindType(int param, std::vector<std::uint8_t> const *blob)
        : param_num_(param), param_str_(nullptr), ptr_(blob),
          type_(SQLITE_BLOB){};

    BindType(int param, std::string const *text)
        : param_num_(param), param_str_(nullptr), ptr_(text),
          type_(SQLITE_TEXT){};

    BindType(const char *param, const std::vector<std::uint8_t> *blob)
        : param_num_(0), param_str_(param), ptr_(blob), type_(SQLITE_BLOB){};

    BindType(const char *param, std::string const *text)
        : param_num_(0), param_str_(param), ptr_(text), type_(SQLITE_TEXT){};
  };


private:
  Database const &db_;
  sqlite3_stmt *stmt_{};// statement
  char *zErrMsg_{};
  int rc_{};
  std::vector<BindType> bind_;
  mutable Row _row = {};
  mutable bool _first_row = false;
};

// wrapper functions ///////////////////////////////////////////////////////////
void create_table(Connection &, const TableDetails &);

void insert_many(Connection &, const String &, const std::vector<Strings> &);

}// namespace util::sqlite

#endif//RNAMAKE_NEW_SQLITE3_CONNECTION_H
