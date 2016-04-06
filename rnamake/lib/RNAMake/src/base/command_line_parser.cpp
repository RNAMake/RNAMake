//
//  command_line_parser.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/1/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "base/command_line_parser.hpp"


void
CommandLineParser::assign_options(
    CommandLineOptions const & cl_options,
    Options & options,
    String prefix) {
    
    auto name = String("");
    for(auto const & cl_opt: cl_options) {
        name = cl_opt->name();
        
        if(prefix.length() > 0) {
            if(name.find(prefix) == 0) {
                name = name.substr(prefix.size() + 2);
            }
        }
        
        if(! options.has_option(name)) { continue; }
        
        
        if     (cl_opt->type() == OptionType::INT) {
            options.set_value(name, cl_opt->get_int());
        }
        else if(cl_opt->type() == OptionType::FLOAT) {
            options.set_value(name, cl_opt->get_float());
        }
        else if(cl_opt->type() == OptionType::BOOL) {
            options.set_value(name, cl_opt->get_bool());
        }
        else if(cl_opt->type() == OptionType::STRING) {
            options.set_value(name, cl_opt->get_string());
        }
    }
    
    
}