//
// Created by Joseph Yesselman on 12/30/17.
//

#ifndef RNAMAKE_NEW_SQLITE3_CONNECTION_H
#define RNAMAKE_NEW_SQLITE3_CONNECTION_H

#include <base/string.hpp>
#include <util/sqlite/database.hpp>
#include <util/sqlite/field.hpp>
#include <util/sqlite/table_details.hpp>

namespace util::sqlite {

class Connection {
public: // construction ///////////////////////////////////////////////////////
  inline explicit Connection(const String & db_name) {
    _db.open(db_name);
  }

  ~Connection() = default;

  Connection(const Connection & other) {
    _db.open(other._db.get_name());
  }

public: // interface to get rows //////////////////////////////////////////////
  int setup_row_iteration(const String &) const;

  const Row &next() const;

  const Row &get_first_row(const String &) const;

  void abort_iterate_rows() const;

  bool has_more_rows() const;

public: // send sql commands //////////////////////////////////////////////////
  int exec(const String &) const;

  void start_transaction() const;

  void commit_transaction() const;

  void rollback_transaction() const;

public: // get other table infos //////////////////////////////////////////////
  TableDetailsOP get_table_details(const String &) const;

  inline const String &get_database_name() const { return _db.get_name(); }

  inline bool is_database_created() const { return _db.is_created(); }

public: // binding ////////////////////////////////////////////////////////////
  // bind a BLOB or TEXT to query
  // CAUTION: vector and string MUST BE constant until end of query execution of
  // exec()/use()/store()
  void bind(int, const Blob &);

  void bind(int, const String &);

  void bind(const char *, const Blob &);

  void bind(const char *, const String &);

private:
  int _execute(const String &) const;

  void _prepare(const String &) const;

  const Row &_generate_row_from_statement() const;

  int _do_bind() const;

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
  Database _db;
  mutable sqlite3_stmt *_stmt = nullptr; // statement
  mutable char *_zErrMsg = new char[1];
  mutable int _rc = 0;
  mutable std::vector<BindType> _bind;
  mutable Row _row = {};
  mutable bool _first_row = false;
  mutable bool _done = true;
};

// wrapper functions ///////////////////////////////////////////////////////////
void create_table(const Connection &, const TableDetails &);

void insert_many(const Connection &, const String &, const std::vector<Strings>
    &);

} // namespace util::sqlite

#endif // RNAMAKE_NEW_SQLITE3_CONNECTION_H
