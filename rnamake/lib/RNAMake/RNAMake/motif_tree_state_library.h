//
//  motif_tree_state_library.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_state_library__
#define __RNAMake__motif_tree_state_library__

#include <stdio.h>
#include "motif_type.h"
#include "motif_tree_state.h"

class MotifTreeStateLibrary {
public:
    MotifTreeStateLibrary(MotifType const &);
    ~MotifTreeStateLibrary() {}

public: //getters:
    
    
private:
    void
    _load_states_from_file(String const &);

private:
    MotifType mtype_;
    MotifTreeStates motif_tree_states_;
};

#endif /* defined(__RNAMake__motif_tree_state_library__) */
