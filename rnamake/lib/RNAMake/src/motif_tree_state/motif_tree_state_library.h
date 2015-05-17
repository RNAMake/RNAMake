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
#include <map>

//RNAMake Headers
#include "motif/motif_type.h"
#include "motif_tree_state/motif_tree_state.h"

class MotifTreeStateLibrary {
public:
    MotifTreeStateLibrary() {}
    MotifTreeStateLibrary(MotifType const &, int new_s=0);
    MotifTreeStateLibrary(MotifType const &, int, int new_s=0);
    MotifTreeStateLibrary(String const &, int new_s=0);
    ~MotifTreeStateLibrary() {}

public:
    MotifTreeStateOP const &
    get_state(String const &);
    
    MotifTreeStateOP const &
    get_state_no_error(String const &);
    
public: //getters:
    inline
    MotifTreeStateOPs const &
    motif_tree_states() const {
        return motif_tree_states_;
    }
    
    inline
    MotifType const &
    mtype() const { return mtype_; }
    
private:
    void
    _load_states_from_file(String const &, int);

private:
    MotifType mtype_;
    MotifTreeStateOPs motif_tree_states_;
    std::map<String, MotifTreeStateOP> motif_tree_state_dict_;
    MotifTreeStateOP null_;
    int new_s_;
};

typedef std::shared_ptr<MotifTreeStateLibrary> MotifTreeStateLibraryOP;
typedef std::vector<MotifTreeStateLibraryOP>   MotifTreeStateLibraryOPs;

#endif /* defined(__RNAMake__motif_tree_state_library__) */
