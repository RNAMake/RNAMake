//
//  sqlite_library.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sqlite_library__
#define __RNAMake__sqlite_library__

#include <stdio.h>
#include <filesystem>

#include "base/types.h"
#include <base/file_io.h>
#include <base/string.h>
#include <base/log.h>

#include <sqlite_modern/sqlite_modern_cpp.h>

namespace resources {

/*
 * Exception for sqlite library
 */


class SqliteLibraryException : public std::runtime_error {
public:
    /**
     * Standard constructor for SqliteLibraryException
     * @param   message   Error message for sqlite libraries
     */
    SqliteLibraryException(String const & message) :
            std::runtime_error(message) {}
};

class SqliteLibrary {
public:
    SqliteLibrary() {}

protected:
    void
    _setup(
            String const &);

    String
    _get_path(
            String const &);

protected:
    StringStringMap libnames_;
    String name_;
    int max_size_;

};

void
build_sqlite_library(String const& , std::vector<Strings>const & , Strings const& , String const& );

void
sqlite3_escape(String &);

void
sqlite3_escape(Strings &);

}

#endif /* defined(__RNAMake__sqlite_library__) */
