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
    auto cmd_opt = CommandLineOption("test", String("test"), OptionType::STRING, false);
   
    auto cl_opts = CommandLineOptions();
    cl_opts.add_option("test", 5, OptionType::INT, false);
    auto val = cl_opts.get_int("test");
    
    //CommandLineOptions cl_opts;
    //cl_opts.add_option("test", "", OptionType::FLOAT, "5", false);
    
    return 1;
}

/*int
CL_OptionUnittest::test_parse_1() {
    
    CommandLineOptions cl_opts;
    cl_opts.add_option("test", "", OptionType::FLOAT, "5", false);
    cl_opts.add_option("test_2", "", OptionType::STRING, "test_3", false);

    auto cla = CommandLineArgs("-test 7.0");
    auto opts = cl_opts.parse_command_line(cla.argc, cla.argv());

    if(opts.get_float("test") != 7.0f) {
        throw UnittestException("did not parse float option correctly");
    }
    
    if(opts.get_string("test_2") != "test_3") {
        throw UnittestException("did not parse string option correctly");
    }
    
    cla = CommandLineArgs("-test 7.0 -test_2 found");
    opts = cl_opts.parse_command_line(cla.argc, cla.argv());

    if(opts.get_string("test_2") != "found") {
        std::cout << opts.get_string("test_2") << std::endl;
        throw UnittestException("did not parse string option correctly");
    }

    return 1;
    
}

int
CL_OptionUnittest::test_parse_2() {
    CommandLineOptions cl_opts;
    cl_opts.add_option("test2", "", OptionType::FLOAT, "5", false);

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
    
    cl_opts.add_option("test", "", OptionType::FLOAT, "5.0", true);
    Options opts = cl_opts.parse_command_line(3, c_argv);
    
    if(opts.get_float("test2") != 5.0f) { return 0; }

    cl_opts.add_option("test3", "", OptionType::FLOAT, "5.0", true);

    try {
        opts = cl_opts.parse_command_line(3, c_argv);
        std::cout << "did not catch exception" << std::endl;
        exit(EXIT_FAILURE);
    } catch(...) {}
    
    return 1;
    
}

int
CL_OptionUnittest::test_parse_3() {
    CommandLineOptions cl_opts;
    cl_opts.add_option("score", "", OptionType::FLOAT, "5", false);
    cl_opts.add_option("sterics", "", OptionType::BOOL, "0", false);

    auto cla = CommandLineArgs("-score 7 -sterics");
    auto opts = cl_opts.parse_command_line(cla.argc, cla.argv());
    auto val = opts.get_bool("sterics");
    if(val != true) {
        throw UnittestException("did not get correct bool option");
    }
    
    
    return 1;
}

int
CL_OptionUnittest::test_add_by_options() {
    auto cl_opts = CommandLineOptions();
    auto opts = Options("TestOptions");
    opts.add_option("test_2", String("test"), OptionType::STRING);
    
    cl_opts.add_options(opts);
    
    char ** argv = new char*[3];
    for(int i = 0; i < 3; i++) {
        argv[i] = new char[100];
    }
    strcpy(argv[0], "");
    strcpy(argv[1], "-test_2");
    strcpy(argv[2], "test_3");
    char const ** c_argv =  ( const char ** ) argv;

    auto opts_2 = cl_opts.parse_command_line(3, c_argv);
    if(opts_2.get_string("test_2") == "test_2") {
        throw UnittestException("did not get expected option value");
    }
    
    return 1;
}
*/

int
CL_OptionUnittest::run() {
    /*if (test_add_option() == 0)  {  std::cout << "test_add_option failed" << std::endl; }
    if (test_parse_1() == 0)     {  std::cout << "test_parse_1 failed" << std::endl; }
    if (test_parse_2() == 0)     {  std::cout << "test_parse_2 failed" << std::endl; }
     */
    return 1;
}



int
CL_OptionUnittest::run_all() {
    String name = "CL_OptionUnittest";
    typedef int (CL_OptionUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_add_option"    ] = &CL_OptionUnittest::test_add_option;
    //func_map["test_parse_1"       ] = &CL_OptionUnittest::test_parse_1;
    //func_map["test_parse_2"       ] = &CL_OptionUnittest::test_parse_2;
    //func_map["test_parse_3"       ] = &CL_OptionUnittest::test_parse_3;
    //func_map["test_add_by_options"] = &CL_OptionUnittest::test_add_by_options;

    int failed = 0;
    for(auto const & kv : func_map) {
        try {
            int result = (this->*kv.second)();
        }
        catch(std::exception const & e) {
            std::cout << name << "::" << kv.first << " returned ERROR! : " << e.what() << std::endl;
        }
    }
    return failed;
}