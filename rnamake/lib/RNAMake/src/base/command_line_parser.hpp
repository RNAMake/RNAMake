//
//  command_line_parser.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/1/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef command_line_parser_hpp
#define command_line_parser_hpp

#include <stdio.h>
#include "base/string.h"
#include "base/option.h"
#include "base/cl_option.h"

class CommandLineParser {
public:
    CommandLineParser() {}
    
    ~CommandLineParser() {}
    
public:
    
    void
    assign_options(
        CommandLineOptions const & cl_options,
        Options & options,
        String prefix = "");
    
    
};

#endif /* command_line_parser_hpp */
