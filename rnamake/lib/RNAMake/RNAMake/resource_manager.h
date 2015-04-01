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
#include "residue_type_set.h"

class ResourceManager {
public:
    static ResourceManager & getInstance() {
        static ResourceManager instance;
        return instance;
    }
    
    
    inline
    ResidueTypeSet const &
    residue_type_set() {
        return rts_;
    }
    
    
private:
    ResourceManager() {
        rts_ = ResidueTypeSet();
    };
    
    
    ResourceManager(ResourceManager const &);
    void operator= (ResourceManager const &);

private:
    ResidueTypeSet rts_;

};




#endif /* defined(__RNAMake__resource_manager__) */
