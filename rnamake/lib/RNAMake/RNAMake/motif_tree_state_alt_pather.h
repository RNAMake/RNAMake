//
//  motif_tree_state_alt_pather.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_state_alt_pather__
#define __RNAMake__motif_tree_state_alt_pather__

#include <stdio.h>
#include <map>
#include "types.h"
#include "settings.h"
#include "motif_tree_state_tree.h"
#include "motif_tree_state_library.h"

class MotifTreeStateAltPather {
public:
    MotifTreeStateAltPather();
    
    ~MotifTreeStateAltPather() {}
    
public:
    std::vector<MotifTreeStateTree> const &
    get_alt_paths(
        MotifTreeStateTree const &);
    
private:
    std::map<String, Strings> sim_dict_;
    std::vector<MotifTreeStateTree> all_trees_;
    MotifTreeStateLibrary mts_lib_;
};

#endif /* defined(__RNAMake__motif_tree_state_alt_pather__) */
