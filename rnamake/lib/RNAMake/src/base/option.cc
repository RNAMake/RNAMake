//
//  Option.cpp
//  REDESIGNC
//
//  Created by Joseph Yesselman on 11/9/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#include <fstream>

//RNAMake Headers
#include "base/string.h"
#include "base/option.h"

void
Option::value(int const & i) {
    if(type_ == OptionType::INT) {
        v_.set<int>(i);
    }
    if(type_ == OptionType::FLOAT) {
        v_.set<float>((float)i);
    }
    
    if(type_ == OptionType::STRING) {
        throw OptionException("attemped to set string but option is a int");
    }
    
    if(type_ == OptionType::BOOL) {
        throw OptionException("attemped to set bool but option is a int");
    }
}

void
Option::value(float const & i) {
    if(type_ == OptionType::INT) {
        v_.set<int>((int)i);
    }
    if(type_ == OptionType::FLOAT) {
        v_.set<float>(i);
    }
    
    if(type_ == OptionType::STRING) {
        throw OptionException("attemped to set string but option is a int");
    }
    
    if(type_ == OptionType::BOOL) {
        throw OptionException("attemped to set bool but option is a int");
    }
}

void
Option::value(String const & i) {
    if(type_ == OptionType::INT) {
        throw OptionException("attemped to set int but option is a String");
    }
    if(type_ == OptionType::FLOAT) {
        throw OptionException("attemped to set float but option is a String");
    }
    
    if(type_ == OptionType::BOOL) {
        throw OptionException("attemped to set bool but option is a String");
    }
    
    if(type_ == OptionType::STRING) {
        v_.set<String>(i);
    }
}

void
Option::value(bool const & i) {
    if(type_ == OptionType::INT) {
        throw OptionException("attemped to set int but option is a bool");
    }
    if(type_ == OptionType::FLOAT) {
        throw OptionException("attemped to set float but option is a bool");
    }
    
    if(type_ == OptionType::STRING) {
        throw OptionException("attemped to set string but option is a bool");
    }
    
    if(type_ == OptionType::BOOL) {
        v_.set<bool>(i);
    }
}














