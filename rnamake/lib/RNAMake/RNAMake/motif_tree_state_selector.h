//
//  motif_tree_state_selector.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/21/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_state_selector__
#define __RNAMake__motif_tree_state_selector__

#include <stdio.h>
#include <map>
#include "types.h"
#include "motif_tree_state_library.h"

class MotifTreeStateSelector {
public:
    MotifTreeStateSelector(
        MotifTreeStateLibraryOPs const &,
        String mode = "helix_flank");
    
    
private:
    MotifTreeStateLibraryOPs mts_libs_;
    std::map<MotifTreeStateLibraryOP, int> lib_map_;
    StringIntMap clash_lists_;
    
};


#endif /* defined(__RNAMake__motif_tree_state_selector__) */
