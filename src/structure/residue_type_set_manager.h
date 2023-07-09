//
//  resource_manager.h
//  RNAMake
//
//  Created by Joseph Yesselman on 3/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__rtsm_resource_manager__
#define __RNAMake__rtsm_resource_manager__

#include <stdio.h>

//RNAMake Headers
#include "structure/residue_type_set.h"

namespace structure {

class ResidueTypeSetManager {
public:
    static ResidueTypeSetManager & getInstance() {
        static ResidueTypeSetManager instance;
        return instance;
    }

public:

    inline
    ResidueTypeSet const &
    residue_type_set() {
        return rts_;
    }


protected:
    ResidueTypeSetManager() { //Prevent construction
        rts_ = ResidueTypeSet();
        //std::cout << "created\n"; should only see this once!!!
    }

    ResidueTypeSetManager(ResidueTypeSetManager const &); //Prevent construction
    void operator=(ResidueTypeSetManager const &);

private:
    ~ResidueTypeSetManager() {}

private:
    ResidueTypeSet rts_;

};

}



#endif /* defined(__RNAMake__resource_manager__) */
