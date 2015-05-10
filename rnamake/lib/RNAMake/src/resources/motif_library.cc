//
//  motif_library.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "util/settings.h"
#include "resources/motif_library.h"

void
MotifLibrary::load_all(int limit) {
    int i = 0;
    connection_.query("SELECT * from motifs");
    Strings values(2);
    while( true) {
        values = connection_.next();
        if (values[0].length() == 0) {
            break;
        }
        MotifOP m (new Motif ( values[0], rts_));
        mdict_[ values[1] ] = m;
        i++;
        if (i > limit ) { break; }
    }

}

MotifOP const
MotifLibrary::get_motif(String const & name) {
    if (mdict_.find(name) != mdict_.end()) {
        return  MotifOP( new Motif(mdict_[name]->copy()));
    }
    
    connection_.query("SELECT * from motifs WHERE name= \'"+name+"\' LIMIT 1");
    Strings values = connection_.next();
    
    if( values[0].length() < 1) {
        throw "cannot find motif with name "+name;
    }
    
    if(mdict_.find(values[0]) == mdict_.end()) {
        mdict_[ values[1] ] = MotifOP(new Motif( values[0], rts_));
    }
    return MotifOP( new Motif(mdict_[ values[1] ]->copy()));
    //return m;

}

bool
MotifLibrary::contains_motif(String const & name) {
    if (mdict_.find(name) != mdict_.end()) { return true; }
    
    connection_.query("SELECT * from motifs WHERE name= \'"+name+"\' LIMIT 1");
    Strings values = connection_.next();
    
    if( values[0].length() < 1) { return false; }
    else { return true; }
    
    
}

String
MotifLibrary::get_db_path(MotifType const & mtype) {
    String path = resources_path() + "/motif_libraries/";
    if     ( mtype == TWOWAY)    { path += "two_ways"; }
    else if( mtype == HELIX )    { path += "helices";  }
    else if( mtype == TCONTACT ) { path += "tertiary_contacts";  }
    else if( mtype == NWAY)      { path += "junctions"; }
    else if( mtype == HAIRPIN)   { path += "hairpins"; }
    else { throw "could not find mtype"; }
    return path+".db";
}


MotifLibrary
unique_twoway_lib() {
    MotifLibrary mlib(TWOWAY);
    String path = resources_path() + "motifs/two_ways/unique_7.dat";
    String line;
    std::ifstream in;
    in.open(path);
    while ( in.good() ) {
        getline(in, line);
        if(line.length() < 10) { break; }
        Strings spl = split_str_by_delimiter(line, " ");
        mlib.get_motif(spl[0]);
    }
    return mlib;
}


MotifLibrary
ideal_helix_lib() {
    MotifLibrary mlib(HELIX);
    std::stringstream ss;
    for(int i = 1; i < 21; i++) {
        ss << "HELIX.LE." << i;
        mlib.get_motif(ss.str());
        ss.str("");
    }
    return mlib;
}





