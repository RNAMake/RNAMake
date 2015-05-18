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

//RNAMake Headers
#include "base/types.h"
#include "base/variant.h"

/*union Values {
    String asString;
    char asChar;
    unsigned char asUChar;
    short asShort;
    unsigned short asUShort;
    int asInt;
    unsigned int asUInt;
    long asLong;
    unsigned long asULong;
    float asFloat;
    double asDouble;

    Values() { asULong = 0; }
    Values(char in) { asUChar = in; }
    Values(unsigned char in) { asChar = in; }
    Values(short in) { asShort = in; }
    Values(unsigned short in) { asUShort = in; }
    Values(int in) { asInt = in; }
    Values(unsigned int in) { asUInt = in; }
    Values(long in) { asLong = in; }
    Values(unsigned long in) { asULong = in; }
    Values(float in) { asFloat = in; }
    Values(double in) { asDouble = in; }
    Values(String in) { asString = in; }
    ~Values() {}

    operator char() { return asChar; }
    operator unsigned char() { return asUChar; }
    operator short() { return asShort; }
    operator unsigned short() { return asUShort; }
    operator int() { return asInt; }
    operator unsigned int() { return asUInt; }
    operator long() { return asLong; }
    operator unsigned long() { return asULong; }
    operator float() { return asFloat; }
    operator double() { return asDouble; }
};*/

class Option {
public:
    using Values = variant<String, float, int>;
    //Option() {}
    
    Option(
        String const & name,
        float const & value):
    name_(name)
    {
        //v_.set<T>(value);

        v_.set<float>(value);

    }
    
    Option(
        String const & name,
        String const & value):
    name_(name)
    {
        v_.set<String>(value);
        
    }
public:
    
    template<typename T>
    T
    value() {
        return v_.get<T>();
    }
    
private:
    String name_;
    Values v_;
    
};


/*class Option {
public:
    Option() {}
    
    Option(
        String const & name,
        float const & value):
    name_(name),
    f_value_(value),
    s_value_(""),
    i_value_(-99),
    data_type_(0)
    {}
    
    Option(
        String const & name,
        String const & value):
    name_(name),
    s_value_(value),
    f_value_(-99),
    i_value_(-99),
    data_type_(1)
    {}
    
    Option(
        String const & name,
        int const & value):
    name_(name),
    i_value_(value),
    data_type_(2)
    {}
    
    template <typename T>
    T
    value() {
        if(data_type_ == 0) { return f_value_; }
        if(data_type_ == 1) { return s_value_; }
        if(data_type_ == 2) { return i_value_; }
        throw "invalid data_type_";
    }
    
    inline
    value(float value) {
        if(data_type_ == 0) { return f_value_; }
        if(data_type_ == 1) { return s_value_; }
        if(data_type_ == 2) { return i_value_; }
        throw "invalid data_type_";
    }
    
public:
    inline
    String const &
    name() { return name_; }

    
    
private:
    String name_;
    int data_type_;
    float f_value_; String s_value_; int i_value_;
    

};*/




/*class Options {
public:
    Options():
    option_map_(std::map<String, Option>())
    {}
    
    ~Options() {}
    
public:
    
    void
    add_option(
        Option opt) {
        option_map_[opt.name()] = opt;
    }
 
    void
    set_option(
        String const & name,
        float const & value) {
        
        option_map_[name].value = value;
        
    }

private:
    bool locked_;
    std::map<String, Option> option_map_;
    
};*/


#endif /* defined(__REDESIGNC__Option__) */
