//
//  FileIO.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 9/29/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __REDESIGNC__FileIO__
#define __REDESIGNC__FileIO__

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>

//RNAMake Headers
#include "base/types.h"


inline
bool
file_exists (String const & name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

inline
int
is_dir(String const & path) {
    if(!file_exists(path)) { return 0; }
    struct stat st;
    lstat(path.c_str(), &st);
    if(S_ISDIR(st.st_mode)) {
        return 1;
    }
    else {
        return 0;
    }
}


Strings
get_lines_from_file(String);




#endif /* defined(__REDESIGNC__FileIO__) */
