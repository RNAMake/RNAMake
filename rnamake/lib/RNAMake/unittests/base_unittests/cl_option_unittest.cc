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
    
    if(val != 5) {
        throw UnittestException("did not get correct value");
    }
    
    return 1;
}

int
CL_OptionUnittest::test_parse_1() {
    
    auto cl_opts = CommandLineOptions();
    cl_opts.add_option("test", 5, OptionType::FLOAT, false);
    cl_opts.add_option("test_2", String("test_3"), OptionType::STRING, false);

    auto cla = CommandLineArgs("-test 7.0");
    cl_opts.parse_command_line(cla.argc, cla.argv());

    if(cl_opts.get_float("test") != 7.0f) {
        throw UnittestException("did not parse float option correctly");
    }
    
    if(cl_opts.get_string("test_2") != "test_3") {
        throw UnittestException("did not parse string option correctly");
    }
    
    if(cl_opts.is_filled("test") == false) {
        throw UnittestException("did not fill option correctly");
    }
    
    cla = CommandLineArgs("-test 7.0 -test_2 found");
    
    auto cl_opts2 = CommandLineOptions();
    cl_opts2.add_option("test", 5, OptionType::FLOAT, false);
    cl_opts2.add_option("test_2", String("test_3"), OptionType::STRING, false);
    
    cl_opts2.parse_command_line(cla.argc, cla.argv());

    if(cl_opts2.get_string("test_2") != "found" ) {
        throw UnittestException("did not parse string option correctly");
    }

    return 1;
    
}

int
CL_OptionUnittest::test_parse_2() {
    auto cl_opts = CommandLineOptions();
    cl_opts.add_option("test2", 5,  OptionType::FLOAT, false);

    auto cla = CommandLineArgs("-test 7.0");
    
    try {
        cl_opts.parse_command_line(cla.argc, cla.argv());
        throw "failed";
    }
    catch(CommandLineOptionException const & e) {}
    catch(...) {
        throw UnittestException("caught unexpected exception");
    }
    
    cl_opts.add_option("test", 5, OptionType::FLOAT, true);
    cla = CommandLineArgs("-test2 7.0");
    
    try {
        cl_opts.parse_command_line(cla.argc, cla.argv());
        throw "failed";
    }
    catch(CommandLineOptionException const & e) { }
    catch(...) {
        throw UnittestException("caught unexpected exception");
    }
    return 1;
    
}

int
CL_OptionUnittest::test_parse_3() {
    CommandLineOptions cl_opts;
    cl_opts.add_option("score", 5, OptionType::FLOAT, false);
    cl_opts.add_option("sterics", false, OptionType::BOOL, false);

    auto cla = CommandLineArgs("-score 7 -sterics");
    cl_opts.parse_command_line(cla.argc, cla.argv());
    auto val = cl_opts.get_bool("sterics");
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
    auto cla = CommandLineArgs("-test_2 test_3");

    cl_opts.parse_command_line(cla.argc, cla.argv());
    if(cl_opts.get_string("test_2") != "test_3") {
        throw UnittestException("did not get expected option value");
    }
    
    return 1;
}

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
    func_map["test_parse_1"       ] = &CL_OptionUnittest::test_parse_1;
    func_map["test_parse_2"       ] = &CL_OptionUnittest::test_parse_2;
    func_map["test_parse_3"       ] = &CL_OptionUnittest::test_parse_3;
    func_map["test_add_by_options"] = &CL_OptionUnittest::test_add_by_options;

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