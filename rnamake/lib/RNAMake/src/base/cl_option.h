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


class CommandLineOption : public Option {
public:
    template<typename T>
    CommandLineOption(
        String const & name,
        T const & value,
        OptionType const & type,
        bool required):
    Option(name, value, type),
    required_(required),
    filled_(false)
    {}
    
    ~CommandLineOption() {}
    
public:
    inline
    bool
    filled() { return filled_; }
    
    inline
    bool
    required() { return required_; }
    
    
public: //setter
    inline
    void
    filled(bool const & filled) { filled_ = filled; }
    
    
private:
    String long_name_;
    bool required_;
    bool filled_;
};

typedef std::shared_ptr<CommandLineOption> CommandLineOptionOP;

class CommandLineOptions{
public:
    CommandLineOptions():
    options_(std::vector<CommandLineOptionOP>())
    {}

public:
    template<typename T>
    void
    add_option(
        String const & name,
        T const & value,
        OptionType const & type,
        bool required) {
        
        auto opt = std::make_shared<CommandLineOption>(name, value, type, required);
        options_.push_back(opt);
    }
    
    void
    add_options(
        Options &);
    
    Options
    parse_command_line(
        int const,
        char const **);
    
public:
    
    inline
    float
    get_int(String const & name) {
        auto opt = _find_option(name);
        return opt->get_int();
    }
    
    inline
    float
    get_float(String const & name) {
        auto opt = _find_option(name);
        return opt->get_float();
    }
    
    inline
    String
    get_string(String const & name) {
        auto opt = _find_option(name);
        return opt->get_string();
    }
    
    inline
    bool
    get_bool(String const & name) {
        auto opt = _find_option(name);
        return opt->get_bool();
    }
    
    inline
    bool
    has_option(
               String const & name) {
        for(auto const & opt : options_) {
            if(opt->name() == name) { return true; }
        }
        
        return false;
    }
    
    template<typename T>
    void
    set_value(
              String const & name,
              T const & val) {
        
        auto opt = _find_option(name);
        opt->value(val);
        
    }
    
    
private:
    CommandLineOptionOP const &
    _find_option(
        String const & name) {
        
        for(auto const & opt : options_) {
            if(opt->name() == name) { return opt; }
        }
        
        throw OptionException("cannot find option with name " + name);
        
    }

    void
    _set_option(
        CommandLineOptionOP const &,
        String const &);
    
private:
    std::vector<CommandLineOptionOP> options_;
   
};


#endif /* defined(__RNAMake__cl_option__) */
