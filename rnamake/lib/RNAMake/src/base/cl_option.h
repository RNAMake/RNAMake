//
//  cl_option.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__cl_option__
#define __RNAMake__cl_option__

#include <stdio.h>
#include <map>
#include <memory>
#include <exception>

//RNAMake Headers
#include "base/types.h"
#include "base/option.h"


struct CommandLineOption {
    
    CommandLineOption(
        String const & ns_name,
        String nl_name = "",
        OptionType notype = OptionType::STRING,
        String nvalue = "",
        bool nrequired = false):
    s_name(ns_name),
    l_name(nl_name),
    otype(notype),
    value(nvalue),
    required(nrequired),
    filled(false)
    {}
    
    String s_name;
    String l_name;
    OptionType otype;
    String value;
    bool required;
    bool filled;
};

typedef std::shared_ptr<CommandLineOption> CommandLineOptionOP;

class CommandLineOptions {
public:
    CommandLineOptions():
    s_cl_opts_(std::map<String, CommandLineOptionOP>()),
    l_cl_opts_(std::map<String, CommandLineOptionOP>())
    {}
    

public:
    void
    add_option(
        String const & s_name,
        String l_name = "",
        OptionType otype = OptionType::STRING,
        String nvalue = "",
        bool required = false);
    
    void
    add_options(
        Options &);
    
    Options
    parse_command_line(
        int const,
        char const **);
    
private:
    Option
    _generate_option(
        CommandLineOptionOP const &,
        String const &);
    
private:
    std::map<String, CommandLineOptionOP> s_cl_opts_;
    std::map<String, CommandLineOptionOP> l_cl_opts_;

};



#endif /* defined(__RNAMake__cl_option__) */
