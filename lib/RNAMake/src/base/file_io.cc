//
//  FileIO.cpp
//  REDESIGNC
//
//  Created by Joseph Yesselman on 9/29/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#include "base/file_io.h"

namespace base {

Strings
get_lines_from_file(String const fname) {
    if(!file_exists(fname)) {
        std::cout << "File: " << fname << "does not exists" << std::endl;
        throw "file Does not Exist";
    }
    
    String line;
    Strings lines;
    std::ifstream input;
    input.open(fname);
    while ( input.good() ) {
        getline(input, line);
        lines.push_back(line);
        
    }
    input.close();
    
    return lines;
}

}