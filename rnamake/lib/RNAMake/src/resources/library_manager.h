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

//RNAMake Headers
#include "base/string.h"
#include "motif/motif_tree.h"
#include "motif/motif_scorer.h"
#include "motif_tree_state/motif_tree_state_library.h"
#include "resources/motif_library.h"

class LibraryManager {
public:
    static LibraryManager & getInstance() {
        static LibraryManager instance;
        return instance;
    }
    
public:
    
    MotifOP
    get_motif(
        String const &,
        int const & end_index = -1,
        String const & end_name = "");
    
    void
    add_motif(String const & path);
    

protected:
    LibraryManager();
    
    LibraryManager(LibraryManager const &);
    void operator=(LibraryManager const &);
    
private:
    ~LibraryManager() {}
    
    void
    _setup_mts_libs();
    
    MotifOP
    _prep_extra_motif_for_asssembly(
        MotifOP const &,
        int end_index = -1,
        String const & end_name = "");
    
    void
    _build_extra_mts(
        MotifOP const &,
        int);

private:
    std::map<String, MotifLibraryOP> mlibs_;
    std::map<String, MotifTreeStateLibraryOP> mts_libs_;
    std::map<String, MotifOP> extra_motifs_;
    std::map<String, MotifTreeStateOP> extra_mts_;
    MotifTree mt_, mt2_;
    MotifScorer ms_;
    
};

#endif /* defined(__RNAMake__library_manager__) */
