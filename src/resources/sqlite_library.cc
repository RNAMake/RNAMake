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

    try {

        auto db = sqlite::database(path);
        auto cmd = std::stringstream {};
        db<<("drop table if exists data_table;\n");
        db<<("create table if not exists data_table("+base::join_by_delimiter(keys," text,")+" primary key (" + primary_key +"));\n");

        cmd<<("insert into data_table VALUES\n");
        for(int ii = 0; ii< data.size(); ++ii)  {
            const auto& entry = data[ii];
            auto line = base::join_by_delimiter(entry,"\',\'");
            line.pop_back();
            line.pop_back();

            cmd<<("(\'"+line+")");

            if(ii < data.size() - 1) {
                cmd<<",";
            } else {
                cmd<<";";
            }
            cmd<<"\n";
        }

        db<<(cmd.str());

    } catch (std::runtime_error const& error ) {
        LOGE<<"Error: "<<error.what();
    }
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
