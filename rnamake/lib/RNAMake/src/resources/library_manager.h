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
    static LibraryManager & getInstance() {
        static LibraryManager instance;
        return instance;
    }
    
public:
    
    MotifOP
    get_motif(String const &);

protected:
    LibraryManager();
    
    LibraryManager(LibraryManager const &);
    void operator=(LibraryManager const &);
    
private:
    ~LibraryManager() {}

private:
    std::map<String, MotifLibraryOP> mlibs_;
    
};

#endif /* defined(__RNAMake__library_manager__) */
