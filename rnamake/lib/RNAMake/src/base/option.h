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

//RNAMake Headers
#include "base/types.h"
#include "base/variant.h"


class Option {
public:
    using Values = variant<String, float, int, bool>;
    Option() {}
    
    Option(
        String const & name,
        float const & value):
    name_(name) {
        v_.set<float>(value);
        type_name_ = String(typeid(float).name());
    }
    
    Option(
        String const & name,
        String const & value):
    name_(name) {
        v_.set<String>(value);
        type_name_ = String(typeid(String).name());
    }
    
    Option(
        String const & name,
        int const & value):
    name_(name) {
        v_.set<int>(value);
        type_name_ = String(typeid(int).name());
    }
    
public:
    
    template<typename T>
    T
    value() {
        return v_.get<T>();
    }
    
    template<typename T>
    void
    value(T const & v) {
        call_type_ = String(typeid(T).name());
        if(call_type_.compare(type_name_) != 0) {
            throw "wrong call type in option\n";
        }
        
        v_.set<T>(v);
    }
    
public: //getters
    
    inline
    String const &
    name() const { return name_; }
    
    
private:
    String name_;
    Values v_;
    String type_name_;
    String call_type_;
};


class Options {
public:
    Options():
    option_map_(std::map<String, Option>())
    {}
    
    ~Options() {}
    
public:
    
    void
    add_option(
        Option const & opt) {
        
        /*if(option_map_.find(opt.name()) != option_map_.end()) {
            throw "trying to overide option in options. not permitted\n";
        }*/
         
        option_map_[opt.name()] = opt;
    }
    
    bool
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
        
    }
    

private:
    bool locked_;
    std::map<String, Option> option_map_;
    
};


#endif /* defined(__REDESIGNC__Option__) */
