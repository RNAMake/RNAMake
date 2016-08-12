//
//  resource_manager.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "resources/resource_manager.h"

MotifOP
RM::motif(
    String const & name,
    String const & end_id,
    String const & end_name) {
    
    for(auto const & kv : mlibs_) {
        if(kv.second->contains(name, end_id, end_name)) {
            return kv.second->get(name, end_id, end_name);
        }
    }
    
    if(added_motifs_.contains(name, end_id, end_name)) {
        return added_motifs_.get(name, end_id, end_name);
    }
    
    throw ResourceManagerException(
        "cannot find motif in resource manager with search: "
        "name=" + name + " end_id=" + end_id + " end_name= " + end_name);
    
}

MotifStateOP
RM::motif_state(
    String const & name,
    String const & end_id,
    String const & end_name) {

    for(auto const & kv : ms_libs_) {
        if(kv.second->contains(name, end_id, end_name)) {
            return kv.second->get(name, end_id, end_name);
        }
    }
    
    if(added_motifs_.contains(name, end_id, end_name)) {
        return added_motifs_.get(name, end_id, end_name)->get_state();
    }
    
    throw ResourceManagerException(
        "cannot find motif state in resource manager with search: "
        "name=" + name + " end_id=" + end_id + " end_name= " + end_name);
}

MotifStateEnsembleOP
RM::motif_state_ensemble(
    String const & name) {
    
    for(auto const & kv : mse_libs_) {
        if(kv.second->contains(name)) {
            return kv.second->get(name);
        }
    }
    
    throw ResourceManagerException(
        "cannot find motif state ensemble in resource manager with search: "
        "name=" + name);
}

void
RM::add_motif(
    String const & path,
    String name) {
    
    auto m = mf_.motif_from_file(path);
    
    if(name != "") { m->name(name); }
    
    int added = 0;
    MotifOPs motifs;
    std::map<Uuid, String, UuidCompare> end_ids;
    for( int i = 0; i < m->ends().size(); i++) {
        auto m_added = mf_.can_align_motif_to_end(m, i);
        if(m_added == nullptr) {
            std::cout << "RESOURCE MANAGER WARNING: cannot create standardized motif for ";
            std::cout << m->name() << " with end" << m->ends()[i]->name() << std::endl;
            continue;
        }
        m_added = mf_.align_motif_to_common_frame(m_added, i);
        if(m_added == nullptr) {
            std::cout << "RESOURCE MANAGER WARNING: cannot create standardized motif for ";
            std::cout << m->name() << " with end" << m->ends()[i]->name() << std::endl;
            continue; }
        
        motifs.push_back(m_added);
        end_ids[m_added->ends()[0]->uuid()] = m_added->end_ids()[0];
    }
    
    if(motifs.size() == 0) {
        throw ResourceManagerException(
            "attempted to add motif from path " + path + " unforunately it has no viable "
            "basepair ends to be build from ");
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
RM::add_motif(
    MotifOP const & m) {
    
    added_motifs_.add_motif(m);
    
    /*int added = 0;
    MotifOPs motifs;
    std::map<Uuid, String, UuidCompare> end_ids;
    for( int i = 0; i < m->ends().size(); i++) {
        auto m_added = mf_.can_align_motif_to_end(m, i);
        if(m_added == nullptr) {
            std::cout << "RESOURCE MANAGER WARNING: cannot create standardized motif for ";
            std::cout << m->name() << " with end" << m->ends()[i]->name() << std::endl;
            continue;
        }
        m_added = mf_.align_motif_to_common_frame(m_added, i);
        if(m_added == nullptr) {
            std::cout << "RESOURCE MANAGER WARNING: cannot create standardized motif for ";
            std::cout << m->name() << " with end" << m->ends()[i]->name() << std::endl;
            continue; }
        
        motifs.push_back(m_added);
        end_ids[m_added->ends()[0]->uuid()] = m_added->end_ids()[0];
    }
    
    if(motifs.size() == 0) {
        throw ResourceManagerException(
            "attempted to add motif with name " + m->name() + " unforunately it has no viable "
            "basepair ends to be build from ");
    }
    
    for(auto & m : motifs) {
        Strings final_end_ids;
        for(auto const & end : m->ends()) {
            final_end_ids.push_back(end_ids[end->uuid()]);
        }
        m->end_ids(final_end_ids);
        added_motifs_.add_motif(m);
    }*/

    
}














