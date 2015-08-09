//
//  sqlite_library.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "util/settings.h"
#include "resources/sqlite_library.h"


String
SqliteLibrary::_get_path(
    String const & libname) {
    
    name_ = libname;
    if(libnames_.find(libname) == libnames_.end()) {
        throw std::runtime_error("cannot find library type in sqlite_library: " + libname);
    }
    return resources_path()+libnames_[libname];
    
}