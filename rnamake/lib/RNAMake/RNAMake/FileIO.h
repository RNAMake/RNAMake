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

#include "types.h"


Strings
split_str_by_delimiter(
	String,
	String);

String
get_base_dir();

inline
bool
file_exists (const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}



#endif /* defined(__REDESIGNC__FileIO__) */
