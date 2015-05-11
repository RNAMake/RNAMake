//
//  string.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "base/string.h"


Strings
split_str_by_delimiter(
    String s,
    String delimiter) {
    
    std::string token;
    std::vector<std::string> tokens;
    size_t pos = 0;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        tokens.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    
    if(s.length() > 0) { tokens.push_back(s); }
    
    return tokens;
    
}

String
filename(
    String const & path) {
    Strings path_spl = split_str_by_delimiter(path, "/");
    return path_spl.back();
}


String
base_dir(
    String const & path) {
    Strings path_spl = split_str_by_delimiter(path, "/");
    String base_path;
    for(int i = 0; i < path_spl.size()-1; i++) {
        base_path += path_spl[i] + "/";
    }
    return base_path;
    
}