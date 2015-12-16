//
//  option_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "option_unittest.h"
#include "base/option.h"
#include <vector>

namespace unittests {

int
OptionUnittest::test_creation() {

    Option opt("test", 5.0f);
    float s = opt.value<float>();

    if((s - 5.0f) > 0.1) {
        return 0;
    }
    
    Option opt2("test", "test");
    String s1 = opt2.value<String>();
    if(s1.compare("test") != 0) {
        return 0;
    }
    
    Option opt3("test", 2);
    int s2 = opt3.value<int>();
    if(s2 != 2) {
        return 0;
    }
    
    try {
        opt2.value<float>();
        std::cout << "did not catch exception" << std::endl;
        exit(EXIT_FAILURE);
    } catch(...) { }

    try {
        opt3.value<float>();
        std::cout << "did not catch exception" << std::endl;
        exit(EXIT_FAILURE);
    } catch(...) { }
    
    
    return 1;
}

int
OptionUnittest::test_add_option() {
    
    Options opts;
    opts.add_option(Option("test", "test"));
    opts.add_option(Option("test_2", 5.0f));
    
    return 1;
}

int
OptionUnittest::test_option() {
    Options opts;
    opts.add_option(Option("test", "test"));
    opts.add_option(Option("test_2", 5.0f));
    
    float v = opts.option<float>("test_2");
    
    if(v != 5.0f) { return 0; }
    
    opts.option<float>("test_2", 7.0f);
    v = opts.option<float>("test_2");

    if(v != 7.0f) { return 0; }
    
    return 1;
}

int
OptionUnittest::run() {
    if (test_creation() == 0)    {  std::cout << "test_creation failed" << std::endl; }
    if (test_add_option() == 0)  {  std::cout << "test_add_option failed" << std::endl; }
    if (test_option() == 0)      {  std::cout << "test_option failed" << std::endl; }

    return 1;
}

int
OptionUnittest::run_all() {
    String name = "OptionUnittest";
    typedef int (OptionUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_creation"   ] = &OptionUnittest::test_creation;
    func_map["test_add_option" ] = &OptionUnittest::test_add_option;
    func_map["test_option"     ] = &OptionUnittest::test_option;
    
    int failed = 0;
    for(auto const & kv : func_map) {
        try {
            int result = (this->*kv.second)();
        }
        catch(...) {
            std::cout << name << "::" << kv.first << " returned ERROR!" << std::endl;
            failed += 1;
        }
    }
    return failed;
}

}