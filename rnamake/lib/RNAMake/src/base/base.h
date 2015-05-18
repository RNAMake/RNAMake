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
    
    inline
    void
    option(
        String const & name,
        float const & value) {
        
        //options_.set_option(name, value);

        
    }
    
protected:
    virtual
    void
    setup_options() = 0;
    
protected:
    //Options options_;
    
};

#endif /* defined(__RNAMake__base__) */
