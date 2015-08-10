//
//  added_motif_library.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/9/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__added_motif_library__
#define __RNAMake__added_motif_library__

#include <stdio.h>

#include "motif/motif.h"
#include "resources/motif_sqlite_library.h"

class AddedMotifLibrary {
public:
    AddedMotifLibrary():
    motifs_(MotifOPs())
    {}
    
    ~AddedMotifLibrary() {}
    
public:
    
    void
    add_motif(
        MotifOP const & m) {
        motifs_.push_back(m);
    }
    
    MotifOP
    get(
        String const & name = dummy_name,
        String const & end_id = dummy_end_id,
        String const & end_name = dummy_name);
    
    MotifOPs
    get_multi(
        String const & name = dummy_name,
        String const & end_id = dummy_end_id,
        String const & end_name = dummy_name);
    
    int
    contains(
        String const & name = dummy_name,
        String const & end_id = dummy_end_id,
        String const & end_name = dummy_name);
    
private:
    
    MotifOPs
    _find_motifs(
        String const & name = dummy_name,
        String const & end_id = dummy_end_id,
        String const & end_name = dummy_name);
    
private:
    MotifOPs motifs_;
};

#endif /* defined(__RNAMake__added_motif_library__) */
