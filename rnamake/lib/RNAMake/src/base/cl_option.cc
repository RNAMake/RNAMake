//
//  cl_option.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "base/cl_option.h"
#include "base/option.h"


void
CommandLineOptions::add_options(
    Options const & opts) {
    
    for(auto const & opt : opts) {
        auto cl_opt = std::make_shared<CommandLineOption>(*opt);
        options_.push_back(cl_opt);
    }
}
 
void
CommandLineOptions::_set_option(
    CommandLineOptionOP const & cl_opt,
    String const & value) {
    
    if(cl_opt->filled()) {
        CommandLineOptionException("commandline option has already has been set!: " +
                                   cl_opt->name());
    }
    
    cl_opt->filled(true);
    
    
    if     (cl_opt->type() == OptionType::STRING) {
        cl_opt->value(String(value));
    }
    else if(cl_opt->type() == OptionType::FLOAT ) {
        cl_opt->value(std::stof(value));
    }
    else if(cl_opt->type() == OptionType::INT   ) {
        cl_opt->value(std::stoi(value));
    }
    else if(cl_opt->type() == OptionType::BOOL  ) {
        cl_opt->value((bool)std::stoi(value));

    }

}

void
CommandLineOptions::parse_command_line(
    int const argc,
    char const ** argv) {
    
    String key = "";
    CommandLineOptionOP cl_opt;
    
    for(int i = 1; i < argc; i++) {
        if(argv[i][0] != '-') { continue; }
        
        key = String(argv[i]);
        key = key.substr(1);

        cl_opt = _find_option(key);
    
    
        if(argc != i + 1 && argv[i+1][0] != '-') {
            _set_option(cl_opt, String(argv[i+1]));
        }
        else {
            _set_option(cl_opt, "1");
        }
            
    }
    
    for(auto const & opt : options_) {
        if(! opt->filled() && opt->required()) {
            throw CommandLineOptionException(opt->name() + " is a required option and was not supplied");
        }
    }
    

}













