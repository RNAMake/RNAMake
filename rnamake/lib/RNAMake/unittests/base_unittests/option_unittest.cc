//
//  option_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "option_unittest.h"
#include "base/option.h"
#include "base/variant.h"
#include <vector>

namespace unittests {

int
OptionUnittest::test_creation() {
    
    auto opt = Option("test", 6, OptionType::FLOAT);
    auto val  = opt.get_float();
    auto val2 = opt.get_int();
    
    try {
        auto val3 = opt.get_string();
        throw "fail";
    }
    catch(OptionException e) {}
    catch(...) {
        throw UnittestException("unexpected error");
    }
    
    opt = Option("test2", String("test"), OptionType::STRING);
    auto val3 = opt.get_string();
    try {
        auto val4 = opt.get_float();
        throw "fail";
    }
    catch(OptionException e) {}
    catch(...) {
        throw UnittestException("unexpected error");
    }
    opt.value("test_2");
    
    opt = Option("test2", false, OptionType::BOOL);
    auto val4 = opt.get_bool();
    opt.value(true);
    
    return 1;
}

int
OptionUnittest::test_add_option() {
    
    auto opts = Options("TestOptions");
    opts.add_option("test", String("test"), OptionType::STRING);
    opts.add_option("test_2", 5, OptionType::INT);
    
    auto val = opts.get_int("test_2");
    if(val != 5) {
        throw UnittestException("did not get expected option value");
    }
    
    return 1;
}

int
OptionUnittest::test_option() {
    
    auto opts = Options("TestOptions");
    opts.add_option("test", String("test"), OptionType::STRING);
    opts.add_option("test_2", 5, OptionType::INT);

    opts.set_value("test_2", 6);
    opts.set_value("test", String("test_2"));

    if(opts.get_string("test") != "test_2") {
        throw UnittestException("did not get expected option value");
    }
    
    if(opts.get_float("test_2") != 6) {
        throw UnittestException("did not get expected option value");
    }
    
    
    
    return 1;
}

int
OptionUnittest::test_iteration() {
    auto opts = Options("TestOptions");
    opts.add_option("test", String("test"), OptionType::STRING);
    opts.add_option("test_2", 5, OptionType::INT);
    
    int count = 0;
    for(auto const & opt : opts) {
        //std::cout << opt.name() << std::endl;
    }
    return 0;
}

int
OptionUnittest::run() {
    
    if (test_creation() == 0)    {  std::cout << "test_creation failed" << std::endl; }
    if (test_add_option() == 0)  {  std::cout << "test_add_option failed" << std::endl; }
    if (test_option() == 0)      {  std::cout << "test_option failed" << std::endl; }
    test_iteration();
    return 1;
}

int
OptionUnittest::run_all() {
    String name = "OptionUnittest";
    typedef int (OptionUnittest::*fptr)();
    std::map<String, fptr> func_map;
    //func_map["test_creation"   ] = &OptionUnittest::test_creation;
    //func_map["test_add_option" ] = &OptionUnittest::test_add_option;
    func_map["test_option"     ] = &OptionUnittest::test_option;
    //func_map["test_iteration"  ] = &OptionUnittest::test_iteration;

    int failed = 0;
    for(auto const & kv : func_map) {
        try {
            int result = (this->*kv.second)();
        }
        catch(std::exception const & e) {
            std::cout << name << "::" << kv.first << " returned ERROR! : " << e.what() << std::endl;
            failed += 1;
        }
    }
    return failed;
}

}