//
//  resource_manager.h
//  RNAMake
//
//  Created by Joseph Yesselman on 3/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__resource_manager__
#define __RNAMake__resource_manager__

#include <stdio.h>

//RNAMake Headers
#include "structure/residue_type_set.h"

class ResourceManager {
public:
    static ResourceManager & getInstance() {
        static ResourceManager instance;
        return instance;
    }
   
public:
    
    inline
    ResidueTypeSet const &
    residue_type_set() {
        return rts_;
    }
    
    
protected:
    ResourceManager() { //Prevent construction
        rts_ = ResidueTypeSet();
        //std::cout << "created\n"; should only see this once!!!
    }
    
    ResourceManager(ResourceManager const &); //Prevent construction
    void operator= (ResourceManager const &);
    
private:
    ~ResourceManager() {}

private:
    ResidueTypeSet rts_;

};




#endif /* defined(__RNAMake__resource_manager__) */
