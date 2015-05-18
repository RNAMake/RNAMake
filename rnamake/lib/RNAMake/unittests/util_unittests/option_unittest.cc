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

int
OptionUnittest::test_creation() {

    Option opt("test", 5.0);
    float s = opt.value<float>();
    std::cout << s << std::endl;
    
    
    Option opt2("test", "test");

    String s1 = opt2.value<String>();
    std::cout << s1 << std::endl;
    /*Options opts;
    opts.add_option(f);
    */
    return 0;
}


int
OptionUnittest::run() {
    test_creation();
    return 1;
}