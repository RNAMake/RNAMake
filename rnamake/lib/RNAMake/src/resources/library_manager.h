//
//  library_manager.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__library_manager__
#define __RNAMake__library_manager__

#include <stdio.h>
#include <map>

#include "base/string.h"
#include "resources/motif_library.h"

class LibraryManager {
public:
    LibraryManager();
    
    ~LibraryManager() {}
    
public:
    
    MotifOP
    get_motif(String const &);
    
private:
    
    std::map<String, MotifLibraryOP> mlibs_;
    
    
};

#endif /* defined(__RNAMake__library_manager__) */
