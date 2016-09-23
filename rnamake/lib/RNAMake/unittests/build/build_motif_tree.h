//
//  build_motif_tree.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__build_motif_tree__
#define __RNAMake__build_motif_tree__

#include <stdio.h>

//RNAMake Headers
#include "motif_data_structures/motif_tree.h"
#include "resources/motif_sqlite_library.h"

class BuildMotifTree {
public:
    BuildMotifTree(
        Strings lib_names = Strings{"ideal_helices", "twoway"}) {
        
        libs_ = std::vector<MotifSqliteLibraryOP>();
        for(auto const & n : lib_names) {
            libs_.push_back(std::make_shared<MotifSqliteLibrary>(n));
        }
    
    }
    
    ~BuildMotifTree() {}
    
public:
    
    MotifTreeOP
    build(
        int size=2);
    
    
    std::unique_ptr<MotifTree>
    build_2(
          int size=2);
private:
    std::vector<MotifSqliteLibraryOP> libs_;
    
};

#endif /* defined(__RNAMake__build_motif_tree__) */
