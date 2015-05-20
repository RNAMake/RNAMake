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
CL_Options::add_option(
    String const & s_name,
    String l_name,
    OptionType otype,
    String nvalue,
    bool required) {
    
    if(s_name.length() == 0 && l_name.length() == 0) {
        throw "cannot add option without a short name or a long name\n";
    }

    CL_OptionOP opt(new CL_Option(s_name, l_name, otype, nvalue, required));
    
    if(s_cl_opts_.find(s_name) != s_cl_opts_.end()) {
        throw "cannot add option " + s_name + " already exists in CL_Options";
    }
    
    s_cl_opts_[s_name] = opt;
    
    if(l_name.length() > 0) {
        l_cl_opts_[l_name] = opt;
    }
    
    
}

Option
CL_Options::_generate_option(
    CL_OptionOP const & cl_opt,
    String const & value) {
    
    cl_opt->filled = true;
    
    String name;
    if(cl_opt->s_name.length() != 0) { name = cl_opt->s_name; }
    else                             { name = cl_opt->l_name; }
        
    if     (cl_opt->otype == STRING_TYPE) { return Option(name, value); }
    else if(cl_opt->otype == FLOAT_TYPE ) { return Option(name, std::stof(value)); }
    else if(cl_opt->otype == INT_TYPE   ) { return Option(name, std::stoi(value)); }
    else { throw "could not generate option from command line"; }

}

Options
CL_Options::parse_command_line(
    int const argc,
    char const ** argv) {
    
    Options opts;
    String key = "";
    CL_OptionOP cl_opt;
    
    int last_arg = argc-1;
    for(int i = 1; i < last_arg; i++) {
        if(argv[i][0] == '-') {
            //
            key = String(argv[i]);
            key = key.substr(1);

            if(argv[i][1] == '-') {
                if(l_cl_opts_.find(key) == l_cl_opts_.end()) {
                    throw std::runtime_error("unknown command line argument: " + key);
                }
                cl_opt = l_cl_opts_[key];
            }
            else {
                if(s_cl_opts_.find(key) == s_cl_opts_.end()) {
                    throw std::runtime_error("unknown command line argument: " + key);
                }
                cl_opt = s_cl_opts_[key];
                
            }
        
          
            Option opt;
            if(argv[i+1][0] != '-') {
                opt = _generate_option(cl_opt, String(argv[i+1]));
            }
            else {
                opt = _generate_option(cl_opt, "0");
            }
            
            opts.add_option(opt);
            
        }
    }
    
    for(auto const & kv : s_cl_opts_) {
        if(kv.second->filled != true) {
            if(kv.second->required == false) {
                Option opt = _generate_option(kv.second, kv.second->value);
                opts.add_option(opt);
            }
            else {
                String message = "missing required argument: ";
                if(cl_opt->s_name.length() != 0 ) { message += "-" + cl_opt->s_name; }
                else                              { message += "--" + cl_opt->l_name; }
                throw message;
            }

        }
    }
    
    return opts;
    
}
















