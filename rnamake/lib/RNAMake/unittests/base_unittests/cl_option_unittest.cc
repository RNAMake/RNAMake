//
//  cl_option_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <string>
#include <cstring>

//RNAMake Headers
#include "cl_option_unittest.h"
#include "base/cl_option.h"

int
CL_OptionUnittest::test_add_option() {
    CL_Options cl_opts;
    cl_opts.add_option("test", "", FLOAT_TYPE, "5.0", false);
    return 1;
}

int
CL_OptionUnittest::test_parse_1() {
    
    CL_Options cl_opts;
    cl_opts.add_option("test", "", FLOAT_TYPE, "5.0", false);
    
    char ** argv = new char*[3];
    for(int i = 0; i < 3; i++) {
        argv[i] = new char[100];
    }
    strcpy(argv[0], "");
    strcpy(argv[1], "-test");
    strcpy(argv[2], "7.0");
    char const ** c_argv =  ( const char ** ) argv;
    
    Options opts = cl_opts.parse_command_line(3, c_argv);
    
    if(opts.option<float>("test") != 7.0f) { return 0; }
    
    return 1;
    
}

int
CL_OptionUnittest::test_parse_2() {
    CL_Options cl_opts;
    cl_opts.add_option("test2", "", FLOAT_TYPE, "5.0", false);

    char ** argv = new char*[3];
    for(int i = 0; i < 3; i++) {
        argv[i] = new char[100];
    }
    strcpy(argv[0], "");
    strcpy(argv[1], "-test");
    strcpy(argv[2], "7.0");
    char const ** c_argv =  ( const char ** ) argv;
    
    try {
        Options opts = cl_opts.parse_command_line(3, c_argv);
        std::cout << "did not catch exception" << std::endl;
        exit(EXIT_FAILURE);
    } catch(...) {}
    
    cl_opts.add_option("test", "", FLOAT_TYPE, "5.0", true);
    Options opts = cl_opts.parse_command_line(3, c_argv);
    
    if(opts.option<float>("test2") != 5.0f) { return 0; }

    cl_opts.add_option("test3", "", FLOAT_TYPE, "5.0", true);

    try {
        opts = cl_opts.parse_command_line(3, c_argv);
        std::cout << "did not catch exception" << std::endl;
        exit(EXIT_FAILURE);
    } catch(...) {}
    
    return 1;
    
}

int
CL_OptionUnittest::run() {
    if (test_add_option() == 0)  {  std::cout << "test_add_option failed" << std::endl; }
    if (test_parse_1() == 0)     {  std::cout << "test_parse_1 failed" << std::endl; }
    if (test_parse_2() == 0)     {  std::cout << "test_parse_2 failed" << std::endl; }

    return 1;
}

int
CL_OptionUnittest::run_all() {
    String name = "CL_OptionUnittest";
    typedef int (CL_OptionUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_creation"   ] = &CL_OptionUnittest::test_add_option;
    func_map["test_add_option" ] = &CL_OptionUnittest::test_parse_1;
    func_map["test_option"     ] = &CL_OptionUnittest::test_parse_2;
    
    for(auto const & kv : func_map) {
        try {
            int result = (this->*kv.second)();
        }
        catch(...) {
            std::cout << name << "::" << kv.first << " returned ERROR!" << std::endl;
        }
    }
    return 0;
}