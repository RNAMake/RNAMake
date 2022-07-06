//
// Created by Joseph Yesselman on 1/3/18.
//

#include <resource_management/sqlite_library.h>

namespace resource_management {

SqliteLibrary::SqliteLibrary(
        String const & db_path,
        String const & table_name):
        db_(util::sqlite::Database(db_path)),
        conn_(util::sqlite::Connection(db_)),
        table_details_(util::sqlite::TableDetails(table_name)) {

    if(conn_.is_database_created()) {
        std::remove(db_path.c_str());
        throw SqliteLibraryException(
                "SqliteLibrary expects an already built sqlite3 data with a given table, was given db_path: " +
                db_path + " which does not exist!");
    }

    table_details_ = *(conn_.get_table_details(table_name));
}


size_t
SqliteLibrary::get_num_of_rows() {
    query_string_ = "SELECT COUNT(*) FROM " + table_details_.name();
    auto & row = conn_.get_first_row(query_string_);
    int count = row[0].get_int();
    return (size_t)count;
}

void
SqliteLibrary::_generate_query(
        Strings const & retrieved_columns,
        StringStringMap const & restraint_col_and_vals) const {

    query_string_ = "SELECT ";
    int i = 0;
    for(auto const & col_name : retrieved_columns) {
        i++;
        if(! _is_valid_key(col_name) ) {
            throw SqliteLibraryException(
                    "attempting to use colunn key: " + col_name + " does not exist in sqlite library");
        }
        query_string_ += col_name + " ";
        if(i != retrieved_columns.size()) { query_string_ += ","; }
    }
    query_string_ += " FROM " + table_details_.name() + " WHERE ";
    i = 0;
    for(auto const & kv : restraint_col_and_vals) {
        i++;
        if(! _is_valid_key(kv.first) ) {
            throw SqliteLibraryException(
                    "attempting to use colunn key: " + kv.first + " does not exist in sqlite library");
        }
        query_string_ += kv.first + "='" + kv.second + "' ";
        if(i != restraint_col_and_vals.size()) {
            query_string_ += "AND ";
        }
    }
}

bool
SqliteLibrary::_is_valid_key(
        String const & key) const {
    for(auto const & col : table_details_) {
        if(col.name == key) { return true; }
    }
    return false;
}

}