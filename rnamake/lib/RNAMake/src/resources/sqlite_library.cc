//
//  sqlite_library.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdexcept>

//RNAMake Headers
#include "base/settings.h"
#include "resources/sqlite_library.h"


String
SqliteLibrary::_get_path(
    String const & libname) {
    
    name_ = libname;
    if(libnames_.find(libname) == libnames_.end()) {
        auto options = String("");
        for(auto const & kv : libnames_) {
            options += kv.first + " ";
        }
        
        throw SqliteLibraryException(
            "cannot find library type in sqlite_library: " + libname + 
            " valid options are: " + options );
    }
    return resources_path()+libnames_[libname];
    
}