//
// Created by Joseph Yesselman on 12/30/17.
//

#ifndef RNAMAKE_NEW_SQLITE_FIELD_H
#define RNAMAKE_NEW_SQLITE_FIELD_H

#include <sqlite3.h>

#include <base/string.hpp>
#include <base/types.hpp>
#include <util/exception.hpp>
#include <utility>

namespace util::sqlite {

typedef std::vector<std::uint8_t> Blob;

class Field {
public:
  inline Field() = default;
  inline Field(String &&name, int i)
      : _name(std::move(name)), _int(i), _type(SQLITE_INTEGER) {}
  inline Field(String &&name, sqlite3_int64 i64)
      : _name(std::move(name)), _int(i64), _type(SQLITE_INTEGER) {}
  inline Field(String &&name, double d)
      : _name(std::move(name)), _float(d), _type(SQLITE_FLOAT) {}
  inline Field(String &&name, String &&s)
      : _name(std::move(name)), _str(std::move(s)), _type(SQLITE_TEXT) {}
  inline Field(String &&name, Blob b)
      : _name(std::move(name)), _blob(std::move(b)), _type(SQLITE_BLOB) {}

  inline ~Field() = default;

public:// getters
  [[nodiscard]] inline int get_int() const {
    if (_type != SQLITE_INTEGER) {
      String msg = "cannot get int value, field is not int type value";
      base::log_and_throw<SqliteException>(msg);
    }
    return static_cast<int>(_int);
  }

  [[nodiscard]] inline sqlite3_int64 get_int64() const {
    if (_type != SQLITE_INTEGER) {
      String msg = "cannot get int value, field is not int type value";
      base::log_and_throw<SqliteException>(msg);
    }
    return _int;
  }

  [[nodiscard]] inline double get_double() const {
    if (_type != SQLITE_FLOAT) {
      String msg = "cannot get double value, field is not double type value";
      base::log_and_throw<SqliteException>(msg);
    }
    return _float;
  }

  [[nodiscard]] inline Blob get_blob() const {
    if (_type != SQLITE_BLOB) {
      String msg = "cannot get blob value, field is not blob type value";
      base::log_and_throw<SqliteException>(msg);
    }
    return _blob;
  }

  [[nodiscard]] inline String get_str() const {
    if (_type != SQLITE_TEXT) {
      String msg = "cannot get str value, field is not str type value";
      base::log_and_throw<SqliteException>(msg);
    }
    return _str;
  }


public:
  [[nodiscard]] inline int get_type() const { return _type; }

  [[nodiscard]] inline String get_name() const { return _name; }

  [[nodiscard]] inline String get_type_str() const {
    if (_type == SQLITE3_TEXT) {
      return _str;
    } else if (_type == SQLITE_INTEGER) {
      return std::to_string(_int);
    } else if (_type == SQLITE_FLOAT) {
      return std::to_string(_float);
    } else if (_type == SQLITE_BLOB) {
      return "BLOB";
    } else {
      throw SqliteException("not a valid sqlite type");
    };
  }


private:
  std::int64_t _int = 0;// int data
  double _float = 0.0f; // float data
  Blob _blob = {};      // vector (blob) data
  String _str;          // string (text) data
  String _name;         // field (col) name
  int _type = 0;        // sqlte type
};

typedef std::vector<Field> Row;
typedef std::shared_ptr<Row> RowOP;

}// namespace util::sqlite


#endif//RNAMAKE_NEW_SQLITE_FIELD_H
