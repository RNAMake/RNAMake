//
// Created by Joseph Yesselman on 1/3/18.
//

#ifndef RNAMAKE_NEW_SQLITE_LIBRARY_H
#define RNAMAKE_NEW_SQLITE_LIBRARY_H

#include <util/sqlite/connection.hpp>

namespace resource_management {


class SqliteLibraryException : public std::runtime_error {
public:
    /**
     * Standard constructor for SqliteLibraryException
     * @param   message   Error message for SqliteLibrary
     */
    explicit SqliteLibraryException(String const & message) :
            std::runtime_error(message) {}
};


class SqliteLibrary {
public:
    SqliteLibrary(
            String const &,
            String const &);

    virtual
    ~SqliteLibrary() = default;

public:
    size_t
    get_num_of_rows();


protected:
    void
    _generate_query(
            Strings const &,
            StringStringMap const &) const;

    bool
    _is_valid_key(
            String const &) const;

protected:
    util::sqlite::Database db_;
    util::sqlite::TableDetails table_details_;
    mutable util::sqlite::Connection conn_; // sqlite3 api commands cannot be const
    mutable String query_string_;
};

}


#endif //RNAMAKE_NEW_SQLITE_LIBRARY_H
