//
//  path_follower.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/18/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__path_follower__
#define __RNAMake__path_follower__

#include <stdio.h>


#include "motif_data_structures/motif_graph.h"

class PathFollower {
public:
    PathFollower() {}
    
    ~PathFollower() {}
        
public:
    
    void
    setup(
        Points const & path);
    
        
private:
    MotifGraphOP mg_;
    Points path_;
        
};

#endif /* defined(__RNAMake__path_follower__) */
