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
#include "types.h"
#include "motif.h"
#include "motif_type.h"
#include "sqlite3_connection.h"
#include "residue_type_set.h"

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
    load_all(int);
    
public:
    
    size_t const
    size() { return mdict_.size(); } 
    
    
private:
    Sqlite3Connection connection_;
    MotifType mtype_;
    StringMotifMap mdict_;
    ResidueTypeSet rts_;
};

#endif /* defined(__RNAMake__motif_library__) */
