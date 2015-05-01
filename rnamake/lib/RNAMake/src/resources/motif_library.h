//
//  motif_library.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_library__
#define __RNAMake__motif_library__

#include <stdio.h>
#include <map>

//RNAMAke Headers
#include "base/types.h"
#include "structure/residue_type_set.h"
#include "motif/motif.h"
#include "motif/motif_type.h"
#include "resources/sqlite3_connection.h"

typedef std::map<String, MotifOP> StringMotifMap;

class MotifLibrary {
public:
    
    MotifLibrary(MotifType const & mtype):
        connection_ ( Sqlite3Connection(mtype) ),
        mtype_ ( mtype ) {
        mdict_ = StringMotifMap();
        rts_ = ResidueTypeSet();
    }
    
    ~MotifLibrary() {}

public:
    
    MotifOP const
    get_motif(String const &);
    
    void
    load_all(int limit = 9999);
    
    bool
    contains_motif(String const &);
    
public:
    
    size_t const
    size() { return mdict_.size(); } 
    
    MotifOPs
    motifs() {
        MotifOPs motifs;
        for(auto const & kv : mdict_) { motifs.push_back(kv.second); }
        return motifs;
    }
    
private:
    Sqlite3Connection connection_;
    MotifType mtype_;
    StringMotifMap mdict_;
    ResidueTypeSet rts_;
};

typedef std::shared_ptr<MotifLibrary> MotifLibraryOP;
typedef std::vector<MotifLibraryOP> MotifLibraryOPs;

MotifLibrary
unique_twoway_lib();

MotifLibrary
ideal_helix_lib();


#endif /* defined(__RNAMake__motif_library__) */
