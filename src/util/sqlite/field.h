//
// Created by Joseph Yesselman on 12/30/17.
//

#ifndef RNAMAKE_NEW_SQLITE_FIELD_H
#define RNAMAKE_NEW_SQLITE_FIELD_H

#include <sqlite3.h>

#include <base/exception.hpp>
#include <base/string.hpp>
#include <base/types.hpp>
#include <utility>

namespace util::sqlite {

class SqliteException : public std::runtime_error {
public:
  explicit SqliteException(const String &message)
      : std::runtime_error(message) {}
};

class Field {
public:
  inline Field() = default;
  inline Field(const String &name, int i)
      : _name(name), _int(i), _type(SQLITE_INTEGER) {}
  inline Field(const String &name, sqlite3_int64 i64)
      : _name(name), _int(i64), _type(SQLITE_INTEGER) {}
  inline Field(const String &name, double d)
      : _name(name), _float(d), _type(SQLITE_FLOAT) {}
  inline Field(const String &name, const String &s)
      : _name(name), _str(s), _type(SQLITE_TEXT) {}
  inline Field(const String &name, const std::vector<std::uint8_t> &v)
      : _name(name), _vec(v), _type(SQLITE_BLOB) {}

public:// getters
  [[nodiscard]] inline int get_int() const {
    if (_type != SQLITE_INTEGER) {
      String msg = "cannot get int value, field is not int type value";
      base::log_and_throw<SqliteException>(msg);
    }
    return static_cast<int>(_int);
  }

  inline explicit operator sqlite3_int64() const {
    if (_type != SQLITE_INTEGER) {
      String msg = "cannot get int value, field is not int type value";
      base::log_and_throw<SqliteException>(msg);
    }
    return _int;
  }

  inline explicit operator double() const {
    if (_type != SQLITE_FLOAT) {
      String msg = "cannot get double value, field is not double type value";
      base::log_and_throw<SqliteException>(msg);
    }
    return _float;
  }

  inline explicit operator String() const { return get_str(); }

  explicit operator std::vector<std::uint8_t>() const { return _vec; }


public:
  [[nodiscard]] inline int get_type() const { return _type; }

  [[nodiscard]] inline String get_name() const { return _name; }

  [[nodiscard]] inline String get_str() const {
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
  std::int64_t _int = 0;              // int data
  double _float = 0.0f;               // float data
  std::vector<std::uint8_t> _vec = {};// vector (blob) data
  String _str;                        // string (text) data
  String _name;                       // field (col) name
  int _type;                          // sqlte type
};

typedef std::vector<Field> Row;
typedef std::shared_ptr<Row> RowOP;

}// namespace util::sqlite


#endif//RNAMAKE_NEW_SQLITE_FIELD_H
