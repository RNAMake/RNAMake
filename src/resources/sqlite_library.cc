//
//  sqlite_library.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdexcept>
#include <sqlite3.h>
//RNAMake Headers
#include "base/settings.h"
#include "resources/sqlite_library.h"

namespace resources {

String
SqliteLibrary::_get_path(
        String const & libname) {

    name_ = libname;
    if (libnames_.find(libname) == libnames_.end()) {
        auto options = String("");
        for (auto const & kv : libnames_) {
            options += kv.first + " ";
        }

        throw SqliteLibraryException(
                "cannot find library type in sqlite_library: " + libname +
                " valid options are: " + options);
    }
    return base::resources_path() + libnames_[libname];

}

void
build_sqlite_library(String const& path, std::vector<Strings>const & data, Strings const& keys, String const& primary_key) {
    if(data.empty()) {
        LOGE<<"Error: no data provided. file write to "<<path<<" aborted";
        return;
    }

    if(data[0].size() != keys.size()) {
        throw SqliteLibraryException("length of each row must be the same length as keys");
    }

    if(!base::file_exists(path)) {
        auto outfile = std::ofstream(path);
        outfile<<'\0';
        outfile.close();
    }

    auto sqlite_command = String{"DROP TABLE IF EXISTS data_table; CREATE TABLE data_table("} + base::join_by_delimiter(keys," TEXT,") + " PRIMARY KEY (" + primary_key + "));";

    sqlite3 *connection = nullptr ;

    auto code = sqlite3_open(path.c_str(), &connection);

    if(code || connection == nullptr) {
        std::cout<<sqlite3_errstr(code)<<std::endl;
        return;
    }

    char* error;

    for(auto&& entry : data)  {

        auto line = base::join_by_delimiter(entry,"\',\'");

        line.pop_back();
        line.pop_back();

        sqlite_command += String{"INSERT INTO data_table VALUES(\'" +line + ");"};
    }

    sqlite3_exec(connection,sqlite_command.c_str(), nullptr, nullptr,&error);
    if(error) {
        LOGE<<"Error ecountered: "<<error<<". For file "<<path;
        sqlite3_free(error);
    }

    sqlite3_close(connection);
}

void
sqlite3_escape(String & unescaped_string ) {
    unescaped_string = base::replace_all(unescaped_string, "\'", "\'\'");
}

void
sqlite3_escape(Strings & unescaped_strings) {
    for(auto& u_string : unescaped_strings) {
        sqlite3_escape(u_string);
    }
}


}
