//
//  resource_manager.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "resources/resource_manager.h"

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
    
    if(added_motifs_.contains(name, end_id, end_name)) {
        return added_motifs_.get(name, end_id, end_name);
    }
    
    throw ResourceManagerException("cannot find motif: ");
    
}

MotifStateOP
ResourceManager::get_state(
    String const & name,
    String const & end_id,
    String const & end_name,
    String const & id) {

    for(auto const & kv : ms_libs_) {
        if(kv.second->contains(name, end_id, end_name, id)) {
            return kv.second->get(name, end_id, end_name, id);
        }
    }
    
    if(added_motifs_.contains(name, end_id, end_name)) {
        return added_motifs_.get(name, end_id, end_name)->get_state();
    }
    
    throw ResourceManagerException("cannot find state: ");
}

MotifStateEnsembleOP
ResourceManager::get_motif_state_ensemble(
    String const & name,
    String const & id) {
    
    for(auto const & kv : mse_libs_) {
        if(kv.second->contains(name, id)) {
            return kv.second->get(name, id);
        }
    }
    
    throw ResourceManagerException("cannot find motif_state_ensemble");
    
    
}

void
ResourceManager::add_motif(
    String const & path) {
    
    auto m = mf_.motif_from_file(path);
    MotifOPs motifs;
    std::map<Uuid, String, UuidCompare> end_ids;
    for( int i = 0; i < m->ends().size(); i++) {
        auto m_added = mf_.can_align_motif_to_end(m, i);
        if(m_added == nullptr) { continue; }
        m_added = mf_.align_motif_to_common_frame(m_added, i);
        if(m_added == nullptr) { continue; }
        motifs.push_back(m_added);
        end_ids[m_added->ends()[0]->uuid()] = m_added->end_ids()[0];
    }
    
    for(auto & m : motifs) {
        Strings final_end_ids;
        for(auto const & end : m->ends()) {
            final_end_ids.push_back(end_ids[end->uuid()]);
        }
        m->end_ids(final_end_ids);
        added_motifs_.add_motif(m);
    }
    
}


void
ResourceManager::register_motif(
    MotifOP const & m) {
    added_motifs_.add_motif(m);
}




