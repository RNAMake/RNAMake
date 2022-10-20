//
//  sqlite_library.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdexcept>
//#include <sqlite3.h>
// RNAMake Headers
#include "base/settings.h"
#include "base/types.hpp"
#include "resources/sqlite_library.h"

namespace resources {

String SqliteLibrary::_get_path(String const &libname) {

  _name = libname;
  if (_libnames.find(libname) == _libnames.end()) {
    auto options = String("");
    for (auto const &kv : _libnames) {
      options += kv.first + " ";
    }

    throw SqliteLibraryException(
        "cannot find library type in sqlite_library: " + libname +
        " valid options are: " + options);
  }
  return base::resources_path() + _libnames[libname];
}

void build_sqlite_library(String const &path, std::vector<Strings> const &data,
                          Strings const &keys, String const &primary_key) {

  if (data.empty()) {
    LOGE << "Error: no data provided. file write to " << path << " aborted";
    return;
  }

  if (data[0].size() != keys.size()) {
    throw SqliteLibraryException(
        "length of each row must be the same length as keys");
  }

  if (std::filesystem::exists(path)) {
    std::filesystem::remove(path);
  }

  try {
    auto insert_str = String{"INSERT INTO data_table("};
    auto table_str = String{"CREATE TABLE data_table("};
    // db<<("drop table if exists data_table;\n");
    for (auto ii = 0; ii < keys.size(); ++ii) {
      table_str += keys[ii] + " TEXT,";
      insert_str += keys[ii] + ",";
    }
    insert_str.pop_back();
    table_str += "PRIMARY KEY (" + primary_key + "));";
    insert_str += ") VALUES ";
    for (int ii = 0; ii < data.size(); ++ii) {
      const auto &entry = data[ii];
      auto line = base::string::join(entry, "\',\'");
      line.pop_back();
      line.pop_back();

      insert_str += "(\'" + line + ")";

      if (ii < data.size() - 1) {
        insert_str += ",";
      } else {
        insert_str += ";";
      }
      insert_str += "\n";
    }

    auto db = sqlite::database(path);
    db << table_str;
    db << insert_str;
  } catch (std::runtime_error const &error) {
    LOGE << "Error: " << error.what();
  }
}

void sqlite3_escape(String &unescaped_string) {
  unescaped_string = base::replace_all(unescaped_string, "\'", "\'\'");
}

void sqlite3_escape(Strings &unescaped_strings) {
  for (auto &u_string : unescaped_strings) {
    sqlite3_escape(u_string);
  }
}
} // namespace resources
