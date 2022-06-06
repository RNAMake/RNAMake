//
// Created by Joseph Yesselman on 12/31/17.
//

#ifndef RNAMAKE_NEW_TABLE_DETAILS_H
#define RNAMAKE_NEW_TABLE_DETAILS_H

#include <base/types.hpp>

#include <util/sqlite/field.h>

namespace util {
namespace sqlite {


class TableDetails {
public:
    struct ColumnDetails {
        String name;
        String type;
        bool is_primary;
    };

public:
    inline
    TableDetails(
            String const & name):
            name_(name),
            columns_(std::vector<ColumnDetails>()) {}

    ~TableDetails() {}

public:
    //iterator
    typedef std::vector<ColumnDetails>::const_iterator const_iterator;

    const_iterator begin() const { return columns_.begin(); }
    const_iterator end() const { return columns_.end(); }

public:
    inline
    ColumnDetails const &
    operator [] (
            Index i) const { return columns_[i]; }

public:

    inline
    size_t
    size() const { return columns_.size(); }

    void
    add_column(
            String const & name,
            String const & type,
            bool is_primary = false) {
        if(!_is_valid_sqlite_type(type)) { throw SqliteException("not a valid sqlite3 type: " + type); }
        if(_name_exists(name)) { throw SqliteException("col name already exists in table cannot repeat"); }

        columns_.push_back(ColumnDetails{name, type, is_primary});
    }

    bool
    has_primary_key() const {
        for(auto const & col : columns_) {
            if(col.is_primary) { return true; }
        }
        return false;
    }

public:
    inline
    String const &
    name() const { return name_; }

private:
    bool
    _is_valid_sqlite_type(
            String const & type) {
        if(type == "TEXT" || type == "BLOB" || type == "INTEGER" || type == "INT" || type == "REAL") { return true; }
        else { return false; }
    }

    bool
    _name_exists(
            String const & name) {
        for(auto const & col : columns_) {
            if(col.name == name) { return true; }
        }
        return false;
    }

private:
    String name_;
    std::vector<ColumnDetails> columns_;

};

typedef std::shared_ptr<TableDetails> TableDetailsOP;

}
}


#endif //RNAMAKE_NEW_TABLE_DETAILS_H
