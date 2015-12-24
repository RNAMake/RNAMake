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
        Option const & opt):
    type_(opt.type_),
    v_(opt.v_),
    name_(opt.name_)
    {}
    
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
        bool const & value,
        OptionType const & type):
    name_(name),
    type_(type) {
        
        
        if(type_ != OptionType::BOOL) {
            throw OptionException("calling wrong constructor for Option");
        }
        
        v_.set<bool>(value);

    
    }
    
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
        
        else if(type_ == OptionType::BOOL) {
            v_.set<bool>(bool(value));
        }
    }
    
public:
    
    float
    get_float() {
        if(type_ == OptionType::STRING) {
            throw OptionException("attemped to get float but option is a string");
        }
        
        if(type_ == OptionType::BOOL) {
            throw OptionException("attemped to get float but option is a bool");
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
        
        if(type_ == OptionType::BOOL) {
            throw OptionException("attemped to get int but option is a bool");
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
        if(type_ == OptionType::BOOL) {
            throw OptionException("attemped to get string but option is a bool");
        }
        
        return v_.get<String>();
    }
    
    bool const &
    get_bool() {
        
        if(type_ == OptionType::INT) {
            throw OptionException("attemped to get bool but option is a int");
        }
        if(type_ == OptionType::FLOAT) {
            throw OptionException("attemped to get bool but option is a float");
        }
        if(type_ == OptionType::STRING) {
            throw OptionException("attemped to get bool but option is a string");
        }
        return v_.get<bool>();
        
    }
    
    void
    value(String const &);
    
    void
    value(float const &);
    
    void
    value(int const &);
    
    void
    value(bool const &);
    
    
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

typedef std::shared_ptr<Option> OptionOP;

class Options {
public:
    Options(String const & name):
    name_(name),
    options_(std::vector<OptionOP>()),
    locked_(false)
    {}
    
    Options(Options const & opts):
    name_(opts.name_),
    locked_(opts.locked_),
    options_(std::vector<OptionOP>()) {
        for(auto const & opt : opts.options_) {
            auto new_opt = std::make_shared<Option>(*opt);
            options_.push_back(new_opt);
        }
        
    }
    
    ~Options() {}
    
public:
    typedef std::vector<OptionOP>::iterator iterator;
    typedef std::vector<OptionOP>::const_iterator const_iterator;
    
    iterator begin() { return options_.begin(); }
    iterator end()   { return options_.end(); }
    
    const_iterator begin() const { return options_.begin(); }
    const_iterator end()   const { return options_.end(); }

    
public:
    
    size_t
    size() { return options_.size(); }
    
    inline
    void
    lock_option_adding() { locked_ = true; }
    
    template<typename T>
    void
    add_option(
        String const & name,
        T const & val,
        OptionType const & type) {
        
        if(locked_) {
            throw OptionException("options are locked you cannot add any new ones");
        }
        
        auto opt = std::make_shared<Option>(name, val, type);
        options_.push_back(opt);
        
    }
    
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
    

    template<typename T>
    void
    set_value(
        String const & name,
        T const & val) {
        
        auto opt = _find_option(name);
        opt->value(val);
        
    }
    
private:
    OptionOP const &
    _find_option(
        String const & name) {
        
        for(auto const & opt : options_) {
            if(opt->name() == name) { return opt; }
        }
        
        throw OptionException("cannot find option with name " + name);
        
    }
    

private:
    bool locked_;
    String name_;
    std::vector<OptionOP> options_;
    
};


#endif /* defined(__REDESIGNC__Option__) */
