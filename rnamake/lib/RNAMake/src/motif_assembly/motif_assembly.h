//
//  motif_assembly.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_assembly__
#define __RNAMake__motif_assembly__

#include <stdio.h>

//RNAMake Headers
#include "base/base.h"
#include "motif/motif.h"
#include "motif/motif_tree.h"
#include "motif/motif_tree_node.h"

class MotifAssembly : public Base {
public:
    MotifAssembly():
    mt_ (MotifTree())
    {}
    
    ~MotifAssembly() {}
    
public:
    
    MotifTreeNodeOP
    add_motif(
        MotifOP const &,
        MotifTreeNodeOP parent = nullptr,
        String end_name = "",
        String parent_name = "");
    
    MotifTreeNodeOP
    add_motif(
        String const &,
        MotifTreeNodeOP parent = nullptr,
        String end_name = "",
        String parent_name = "");
    
protected:
    
    void
    setup_options();
    
    
private:
    MotifTree mt_;
    
};

#endif /* defined(__RNAMake__motif_assembly__) */
