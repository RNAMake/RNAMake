//
//  resource_manager.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__resource_manager__
#define __RNAMake__resource_manager__

#include <stdio.h>

//RNAMake Headers
#include "resources/motif_sqlite_library.h"


class ResourceManger {
public:
    static ResourceManger & getInstance() {
        static ResourceManger instance;
        return instance;
    }
    
public:
    
    
protected:
    ResourceManger() { //Prevent construction
        
    
    }
    
    ResourceManger(ResourceManger const &); //Prevent construction
    void operator= (ResourceManger const &);
    
private:
    ~ResourceManger() {}
    
private:

    
};


#endif /* defined(__RNAMake__resource_manager__) */
