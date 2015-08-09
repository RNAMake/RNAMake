//
//  motif_sqlite_library.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_sqlite_library__
#define __RNAMake__motif_sqlite_library__

#include <stdio.h>
#include <iostream>

#include "motif/motif.h"
#include "resources/motif_sqlite_connection.h"
#include "resources/sqlite_library.h"

static String dummy_name_ = "", dummy_end_id_ ="", dummy_end_name_ = "";

class MotifSqliteLibrary : public SqliteLibrary {
public:
    
    MotifSqliteLibrary(
        String const & libname) {
        
        libnames_ = get_libnames();
        auto path = _get_path(libname);
        MotifSqliteConnection conn(path);
        connection_ = conn;
        max_size_ = connection_.count();
        
    }

    ~MotifSqliteLibrary() {}
    
public:
    
    static
    StringStringMap
    get_libnames();
    
    MotifOP
    get(
        String const & name = dummy_name_,
        String const & end_id = dummy_end_id_,
        String const & end_name = dummy_name_);

private:
    
    String
    _generate_query(
        String const &,
        String const &,
        String const &);
    
    
private:

    MotifSqliteConnection connection_;
    std::map<String, MotifOP> data_;
    
    
};


#endif /* defined(__RNAMake__motif_sqlite_library__) */
