//
//  added_motif_library.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/9/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "resources/added_motif_library.h"

MotifOPs
AddedMotifLibrary::_find_motifs(
    String const & name,
    String const & end_id,
    String const & end_name) {
    
    MotifOPs motifs;
    for(auto const & m : motifs_) {
        if(name.length() > 0 && name != m->name()) { continue; }
        if(end_id.length() > 0 && end_id != m->end_ids()[0]) { continue; }
        if(end_name.length() > 0 && end_name != m->ends()[0]->name()) { continue; }
        motifs.push_back(m);
    }
    return motifs;
}

MotifOP
AddedMotifLibrary::get(
    String const & name,
    String const & end_id,
    String const & end_name) {
    
    auto motifs =  _find_motifs(name, end_id, end_name);
    if(motifs.size() == 0) {
        throw std::runtime_error("called get in AddedMotifLibrary but returned no motifs");
    }
    return motifs[0];
    
}

MotifOPs
AddedMotifLibrary::get_multi(
    String const & name,
    String const & end_id,
    String const & end_name) {
    
    return _find_motifs(name, end_id, end_name);
}

int
AddedMotifLibrary::contains(
    String const & name,
    String const & end_id,
    String const & end_name) {
    
    auto motifs = _find_motifs(name, end_id, end_name);
    if (motifs.size() == 0) { return 0; }
    else                    { return 1; }
    
}