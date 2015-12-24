//
//  base.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__base__
#define __RNAMake__base__

#include <stdio.h>

//RNAMake Headers
#include "base/option.h"

class Base {
public:
    
    /*template <typename T>
    inline
    void
    option(
        String const & name,
        T const & value) {
        
        options_.option<T>(name, value);
        
        update_var_options();
        
    }
    
    template <typename T>
    inline
    T
    option(
        String const & name) {

        return options_.option<T>(name);
        
    }
    
    inline
    bool
    contain_option(
        String const & name) {
     
        return options_.contains_option(name);
        
    }
    */
    
    
protected:
    virtual
    void
    setup_options() = 0;
    
    virtual
    void
    update_var_options() { }
    
protected:
    Options options_;
    
};

#endif /* defined(__RNAMake__base__) */
