//
//  build_motif_graph.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__build_motif_graph__
#define __RNAMake__build_motif_graph__

#include <stdio.h>

//RNAMake Headers
#include "resources/motif_sqlite_library.h"
#include "motif_data_structures/motif_graph.h"

namespace unittests {

class BuildMotifGraph {
public:
    BuildMotifGraph(
        Strings lib_names = Strings{"ideal_helices", "twoway"}) {
        
        libs_ = std::vector<MotifSqliteLibraryOP>();
        for(auto const & n : lib_names) {
            libs_.push_back(std::make_shared<MotifSqliteLibrary>(n));
        }
    }
    
    ~BuildMotifGraph() {}
    
public:
    
    MotifGraphOP
    build(int size=2);
    
private:
    std::vector<MotifSqliteLibraryOP> libs_;
    
};

}
#endif /* defined(__RNAMake__build_motif_graph__) */
