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

//RNAMake Headers
#include "base/types.h"
#include "base/option.h"

enum OptionType {
    BOOL_TYPE,
    INT_TYPE,
    STRING_TYPE,
    FLOAT_TYPE
};

struct CL_Option {
    
    CL_Option(
        String const & ns_name,
        String nl_name = "",
        OptionType notype = STRING_TYPE,
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

typedef std::shared_ptr<CL_Option> CL_OptionOP;

class CL_Options {
public:
    CL_Options():
    s_cl_opts_(std::map<String, CL_OptionOP>()),
    l_cl_opts_(std::map<String, CL_OptionOP>())
    {}
    

public:
    void
    add_option(
        String const & s_name,
        String l_name = "",
        OptionType otype = STRING_TYPE,
        String nvalue = "",
        bool required = false);
    
    Options
    parse_command_line(
        int const,
        char const **);
    
private:
    Option
    _generate_option(
        CL_OptionOP const &,
        String const &);
    
private:
    std::map<String, CL_OptionOP> s_cl_opts_;
    std::map<String, CL_OptionOP> l_cl_opts_;

};



#endif /* defined(__RNAMake__cl_option__) */
