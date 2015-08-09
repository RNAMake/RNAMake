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

#include "base/types.h"

class SqliteLibrary {
public:
    SqliteLibrary()
    {}
    
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

#endif /* defined(__RNAMake__sqlite_library__) */
