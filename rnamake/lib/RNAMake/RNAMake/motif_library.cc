//
//  motif_library.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_library.h"

//Motif const
//MotifLibrary::get_motif();

void
MotifLibrary::load_all(int limit=9999) {
    int i = 0;
    connection_.query("SELECT * from motifs");
    Strings values(2);
    while( true) {
        values = connection_.next();
        Motif m ( values[0], rts_);
        mdict_[ values[1] ] = m;
        i++;
        if (i > limit ) { break; }
    }

}

Motif const
MotifLibrary::get_motif(String const & name) {
    if (mdict_.find(name) != mdict_.end()) {
        return mdict_[name].copy();
    }
    
    connection_.query("SELECT * from motifs WHERE name= \'"+name+"\' LIMIT 1");
    Strings values = connection_.next();
    
    if( values[0].length() < 1) {
        throw "cannot find motif with name "+name;
    }
    
    Motif m ( values[0], rts_);
    mdict_[ values[1] ] = m;
    return m.copy();
    //return m;

}







