//
// Created by Joseph Yesselman on 12/30/17.
//

#ifndef RNAMAKE_NEW_SQLITE_FIELD_H
#define RNAMAKE_NEW_SQLITE_FIELD_H

#include <sqlite3.h>

#include <base/assertions.h>
#include <base/types.hpp>
#include <base/string.hpp>
#include <base/vector_container.h>

namespace util {
namespace sqlite {

class SqliteException : public std::runtime_error {
public:
    SqliteException(String const & message) :
            std::runtime_error(message) {}
};

class Field {
public:
    inline
    Field(
            String const & name,
            int i):
            name_(name),
            int_(i),
            type_(SQLITE_INTEGER) {}
    inline
    Field(
            String const & name,
            sqlite3_int64 i64):
            name_(name),
            int_(i64),
            type_(SQLITE_INTEGER) {}
    inline
    Field(
            String const & name,
            double d):
            name_(name),
            float_(d),
            type_(SQLITE_FLOAT) {}

    inline
    Field(
            String const & name,
            String const & s):
            name_(name),
            str_(s),
            type_(SQLITE_TEXT) {}

    inline
    Field(
            String const & name,
            std::vector<std::uint8_t> const & v):
            name_(name),
            vec_(v),
            type_(SQLITE_BLOB) {}

public: // getters
    inline
    operator
    int() const {
        expects<SqliteException>(
                type_ == SQLITE_INTEGER,
                "cannot get int value, field is not int type value");
        return static_cast<int>(int_);
    }

    inline
    operator
    sqlite3_int64() const {
        expects<SqliteException>(
                type_ == SQLITE_INTEGER,
                "cannot get int value, field is not int type value");
        return int_;
    }

    inline
    operator
    double() const {
        expects<SqliteException>(
                type_ == SQLITE_FLOAT,
                "cannot get double value, field is not double type value");
        return float_;
    }

    inline
    operator
    String() const {
        return get_str();
    }

    operator
    std::vector<std::uint8_t>() const { return vec_; }


public:
    inline
    int
    get_type() const { return type_; }

    inline
    String
    get_name() const { return name_; }

    inline
    String
    get_str() const {
        if     (type_ == SQLITE3_TEXT)   { return str_; }
        else if(type_ == SQLITE_INTEGER) { return std::to_string(int_); }
        else if(type_ == SQLITE_FLOAT)   { return std::to_string(float_); }
        else if(type_ == SQLITE_BLOB)    { return "BLOB"; }
        else {
           throw SqliteException("not a valid sqlite type");
        };
    }


private:
    std::int64_t int_;   // int data
    double float_; // float data
    std::vector<std::uint8_t> vec_;   // vector (blob) data
    String str_;   // string (text) data
    String name_;  // field (col) name
    int type_;  // sqlte type
};

typedef base::VectorContainer<Field> Row;
typedef std::shared_ptr<Row>         RowOP;

}
}


#endif //RNAMAKE_NEW_SQLITE_FIELD_H
