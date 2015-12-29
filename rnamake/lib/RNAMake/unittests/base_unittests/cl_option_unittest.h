//
//  cl_option_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__cl_option_unittest__
#define __RNAMake__cl_option_unittest__

#include <stdio.h>

#include "base/string.h"
#include "unittest.h"

struct CommandLineArgs {
    char ** argv_;
    int argc;
    
    CommandLineArgs(
        String const & s) {
        
        Strings spl = split_str_by_delimiter(s, " ");
        argv_ = new char*[spl.size()+1];
        argv_[0] = new char[20];
        strcpy(argv_[0], "program_name");
        for(int i = 0; i < spl.size(); i++) {
            argv_[i+1] = new char[spl[i].size()];
            strcpy(argv_[i+1], spl[i].c_str());
        }
        argc = (int)(spl.size() + 1);

    }
    
    char const **
    argv() {
        return ( const char ** ) argv_;

    }
    
};

class CL_OptionUnittest : public Unittest {
public:
    CL_OptionUnittest() {}
    
    ~CL_OptionUnittest() {}
    
    int
    size() { return 5; }
    
public:
    
    int
    test_add_option();
    
    int
    test_parse_1();
    
    int
    test_parse_2();
    
    int
    test_parse_3();
    
    int
    test_add_by_options();
    

public:
    
    int
    run();
    
    int
    run_all();
    
};

#endif /* defined(__RNAMake__cl_option_unittest__) */
