//
//  resource_manager.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "resource_manager.h"


MotifOP
ResourceManager::get_motif(
    String const & name,
    String const & end_id,
    String const & end_name,
    String const & id) {
    
    for(auto const & kv : mlibs_) {
        if(kv.second->contains(name, end_id, end_name, id)) {
            return kv.second->get(name, end_id, end_name, id);
        }
    }
    
    throw ResourceManagerException("cannot find motif: ");
    
}