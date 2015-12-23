//
//  Option.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 11/9/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __REDESIGNC__Option__
#define __REDESIGNC__Option__

#include <stdio.h>
#include <vector>
#include <map>
#include <typeinfo>
#include <iostream>
#include <memory>
#include <sstream>

//RNAMake Headers
#include "base/types.h"
#include "base/variant.h"

enum class OptionType {
    BOOL,
    INT,
    STRING,
    FLOAT
};


class OptionException : public std::runtime_error {
public:
    OptionException(
        String const & message):
    std::runtime_error(message)
    {}
};

class Option {
public:
    using Values = variant<String, float, int, bool>;
    Option() {}
    
    Option(
        String const & name,
        float const & value,
        OptionType const & type):
    name_(name),
    type_(type) { v_.set<float>(value); }
    
    Option(
        String const & name,
        String const & value,
        OptionType const & type):
    name_(name),
    type_(type) { v_.set<String>(value); }
    
    Option(
        String const & name,
        int const & value,
        OptionType const & type):
    name_(name),
    type_(type) {
        if(type_ == OptionType::INT) {
            v_.set<int>(value);
        }
        else if(type_ == OptionType::FLOAT) {
            v_.set<float>(float(value));
        }
    }
    
public:
    
    float
    get_float() {
        if(type_ == OptionType::STRING) {
            throw OptionException("attemped to get float but option is a string");
        }
        
        if(type_ == OptionType::INT) {
            return v_.get<int>();
        }
        
        return v_.get<float>();
    }
    
    int
    get_int() {
        if(type_ == OptionType::STRING) {
            throw OptionException("attemped to get int but option is a string");
        }
        
        if(type_ == OptionType::FLOAT) {
            return v_.get<float>();
        }
        
        return v_.get<int>();
    }
    
    String const &
    get_string() {
        
        if(type_ == OptionType::INT) {
            throw OptionException("attemped to get string but option is a int");
        }
        if(type_ == OptionType::FLOAT) {
            throw OptionException("attemped to get string but option is a float");
        }
        
        return v_.get<String>();
    }
    
    template<typename T>
    void
    value(T const & v) {
    }
    
    
public: //getters
    
    inline
    String const &
    name() const { return name_; }
    
    inline
    OptionType const &
    type() const { return type_; }
    
private:
    String name_;
    Values v_;
    OptionType type_;
    
};


class Options {
public:
    Options(String const & name):
    name_(name),
    options_(std::vector<Option>())
    {}
    
    ~Options() {}
    
public:
    typedef std::vector<Option>::iterator iterator;
    typedef std::vector<Option>::const_iterator const_iterator;
    
    iterator begin() { return options_.begin(); }
    iterator end()   { return options_.end(); }
    
    const_iterator begin() const { return options_.begin(); }
    const_iterator end()   const { return options_.end(); }

    
public:
    
    size_t
    size() { return options_.size(); }
    
    void
    add_option(
        Option const & opt) {
        
        options_.push_back(opt);
        
    }
    
    inline
    float
    get_int(String const & name) {
        auto opt = _find_option(name);
        return opt.get_int();
    }
    
    inline
    float
    get_float(String const & name) {
        auto opt = _find_option(name);
        return opt.get_float();
    }
    
    inline
    String
    get_string(String const & name) {
        auto opt = _find_option(name);
        return opt.get_string();
    }
    

    template<typename T>
    void
    set_value(
        String const & name,
        T const & val) {
        
        auto opt = _find_option(name);
        opt.value(val);
        
    }
    
    /*bool
    contains_option(
        String const & name) {
        
        if(option_map_.find(name) != option_map_.end()) {
            return true;
        }
        else {
            return false;
        }
    }

    template <typename T>
    T
    option(
        String const & name) {
        
        if(option_map_.find(name) == option_map_.end()) {
            throw "trying to access option that does not exist\n";
        }
        
        return option_map_[name].value<T>();
        
    }
           
           
    template <typename T>
    void
    option(
        String const & name,
        T const & v) {
        
        if(option_map_.find(name) == option_map_.end()) {
            throw "trying to access option that does not exist\n";
        }
        
        option_map_[name].value(v);
        
    }*/
    
private:
    Option const &
    _find_option(
        String const & name) {
        
        for(auto const & opt : options_) {
            if(opt.name() == name) { return opt; }
        }
        
        throw OptionException("cannot find option with name " + name);
        
    }
    

private:
    bool locked_;
    String name_;
    std::vector<Option> options_;
    
};


#endif /* defined(__REDESIGNC__Option__) */
