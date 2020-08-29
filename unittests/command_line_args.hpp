#ifndef __command_line_args
#define __command_line_args

#include "base/string.h"
#include <iostream>

struct CommandLineArgs {
    char ** argv_;
    int argc;
    
    CommandLineArgs(
        String const & s) {
        Strings spl = base::split_str_by_delimiter(s, " ");
        argv_ = new char*[spl.size()+1];
        argv_[0] = new char[20];
        strcpy(argv_[0], "program_name\0");
        
        for(int i = 0; i < spl.size(); i++) {
            argv_[i+1] = new char[spl[i].size() + 1];
            strcpy(argv_[i+1], spl[i].c_str());
            argv_[i+1][spl[i].size()] = '\0';
        }
        
        argc = (int)(spl.size() + 1);
    }
    
    char const **
    argv() {
        return ( const char ** ) argv_;
        
    }
    
};

#endif
