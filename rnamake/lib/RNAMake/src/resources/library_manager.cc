//
//  library_manager.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <sstream>

//RNAMake Headers
#include "util/settings.h"
#include "structure/basepair_state.h"
#include "motif/motif_type.h"
#include "motif/motif_tree.h"
#include "motif/motif.h"
#include "motif_tree_state/motif_tree_state.h"
#include "resources/motif_library.h"
#include "resources/library_manager.h"

void
LibraryManager::_setup_mts_libs() {
    String path = base_dir() + "/rnamake/unittests/";
    String path1 = resources_path() + "prediction/";
    
    mts_libs_["TWOWAY"  ] = MotifTreeStateLibraryOP(new MotifTreeStateLibrary(path+"test_twoway.new.me",1));
    mts_libs_["HELIX"   ] = MotifTreeStateLibraryOP(new MotifTreeStateLibrary(path+"test_helix.new.me",1));
    mts_libs_["BP_STEPS"] = MotifTreeStateLibraryOP(new MotifTreeStateLibrary(path1+"all.new.me",1));
    
}

LibraryManager::LibraryManager() {
    mlibs_ = std::map<String, MotifLibraryOP>();
    mts_libs_ = std::map<String, MotifTreeStateLibraryOP>();
    ms_ = MotifScorer();
    
    extra_motifs_ = std::map<String, MotifOP>();
    extra_mts_    = std::map<String, MotifTreeStateOP>();
    mlibs_["HELIX"   ] = MotifLibraryOP(new MotifLibrary(HELIX));
    mlibs_["NWAY"    ] = MotifLibraryOP(new MotifLibrary(NWAY));
    mlibs_["TWOWAY"  ] = MotifLibraryOP(new MotifLibrary(TWOWAY));
    mlibs_["TCONTACT"] = MotifLibraryOP(new MotifLibrary(TCONTACT));
    
    String path = resources_path() + "motif_libraries/bp_steps.db";
    MotifLibraryOP mlib (new MotifLibrary(path, HELIX));
    String path2 = resources_path() + "motif_libraries/seq_and_ss";
    MotifLibraryOP mlib2 (new MotifLibrary(path, HELIX));

    mlibs_["BP_STEPS"] = mlib;
    mlibs_["SEQ_AND_SS"] = mlib2;
    
    _setup_mts_libs();
    mt_  = MotifTree();
    mt2_ = MotifTree();
    mt_.add_motif(get_motif("HELIX.IDEAL.6"), nullptr, 1, -1, 0);
    mt_.increase_level();
    
}

MotifOP
LibraryManager::get_motif(
    String const & name,
    int const & end_index,
    String const & end_name) {
    for(auto const & kv : mlibs_) {
        if(kv.second->contains_motif(name)) {
            return kv.second->get_motif(name);
        }
    }
    
    if(extra_motifs_.find(name) != extra_motifs_.end()) {
        if(end_index == -1 && end_name.length() == 0) {
            return extra_motifs_[name];
        }
        
        MotifOP m = extra_motifs_[name];
        return _prep_extra_motif_for_asssembly(m);
    }
    
    throw "cannot find motif " + name + "\n";
    
}

void
LibraryManager::add_motif(String const & path) {
    String name = filename(path);
    extra_motifs_[name] = MotifOP(new Motif(path));
    
}

MotifOP
LibraryManager::_prep_extra_motif_for_asssembly(
    MotifOP const & m,
    int end_index,
    String const & end_name) {
    
    if(end_index == -1 && end_name.length() == 0) {
        throw "cannot call _prep_extra_motif_for_asssembly without end_name or end_index";
    }
    
    if(end_name.length() != 0) {
        BasepairOP end = m->get_basepair_by_name(end_name);
        end_index = m->end_index(end);
    }
    
    MotifTreeNodeOP node = mt_.add_motif(m, nullptr, end_index);
    if(node == nullptr) {
        throw "cannot prep motif for assembly, motif: " + m->name();
    }
    
    MotifOP h_m = get_motif("HELIX.IDEAL.3");
    for(auto const & e : node->available_ends()) {
        MotifTreeNodeOP n = mt_.add_motif(h_m, node, 1, -1, 0);
        if(n == nullptr) {
            e->flip();
        }
        else {
            mt_.remove_last_node();
        }
    }
    
    MotifOP m_copy = mt_.nodes()[2]->motif();
    mt_.remove_node_level();
    
    return m;
    
    
}

void
LibraryManager::_build_extra_mts(
    MotifOP const & m,
    int end_index) {
    
    mt2_.add_motif(m, nullptr, end_index, -1, 0);
    MotifOP m_copy = mt2_.nodes()[1]->motif();
    mt2_.remove_node_level();
    
    BasepairStateOPs ends(m_copy->ends().size());
    int i = 0;
    for(auto const & end : m_copy->ends()) {
        if(end_index == i) { continue; }
        ends[i] = BasepairStateOP(new BasepairState(end->state()));
        i++;
    }
    
    Points beads;
    for(auto const & b : m_copy->beads()) {
        if (b.btype() != PHOS) { beads.push_back(b.center()); }
    }
    
    std::stringstream ss;
    ss << m_copy->name() << "-" << end_index;
    String build_str = m_copy->to_str();
    float score = ms_.score(m);
    MotifTreeStateOP mts(new MotifTreeState(ss.str(), end_index, (int)m->residues().size(), score,
                                            beads, ends, 0, m_copy->to_str()));
    
    extra_mts_[ss.str()] = mts;
}


















