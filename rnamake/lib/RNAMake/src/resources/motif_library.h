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
#include "resources/motif_sqlite_connection.h"

typedef std::map<String, MotifOP> StringMotifMap;

class MotifLibrary {
public:
    
    MotifLibrary(MotifType const & mtype):
        connection_ ( MotifSqliteConnection(get_db_path(mtype))),
        mtype_ ( mtype ) {
        mdict_ = StringMotifMap();
        rts_ = ResidueTypeSet();
    }
    
    MotifLibrary(String const & path, MotifType const & mtype = UNKNOWN):
        connection_(MotifSqliteConnection(path)),
        mtype_(mtype) {
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
    
    String
    get_db_path(MotifType const & mtype);
    
private:
    MotifSqliteConnection connection_;
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
