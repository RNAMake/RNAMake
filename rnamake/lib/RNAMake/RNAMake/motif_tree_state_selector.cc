//
//  motif_tree_state_selector.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/21/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_state_selector.h"
#include "motif_tree_state_search_node.h"
#include "motif_tree_state_library.h"
#include "settings.h"
#include "FileIO.h"

MotifTreeStateSelector::MotifTreeStateSelector(
    MotifTreeStateLibraryOPs const & mts_libs,
    String mode):
    nodes_ ( SelectorNodeOPs() ),
    clash_lists_ ( std::map<String, StringIntMapOP>() ),
    mts_type_pairs_ ( MTSTypePairs() )
{
    
    if     (mode.compare("all") == 0) {
        int i = -1;
        for(auto const & mts_lib : mts_libs) {
            i++;
            SelectorNodeOP node (new SelectorNode(mts_lib, i));
            nodes_.push_back(node);
        }
        for (auto const & n1 : nodes_) {
            for (auto const & n2 : nodes_) {
                n1->add_connection(n2);
            }
        }
    }
    
    else if(mode.compare("round_robin") == 0) {
        int i = -1;
        for(auto const & mts_lib : mts_libs) {
            i++;
            SelectorNodeOP node (new SelectorNode(mts_lib, i));
            nodes_.push_back(node);
        }
        for (auto const & n1 : nodes_) {
            for (auto const & n2 : nodes_) {
                if(n1 == n2) { continue; }
                n1->add_connection(n2);
            }
        }
    }
    
    else if(mode.compare("helix_flank") == 0) {
        MotifTreeStateLibraryOP hlib ( new MotifTreeStateLibrary (HELIX));
        SelectorNodeOP node (new SelectorNode(hlib, 0));
        nodes_.push_back(node);
        int i = -1;
        for( auto const & mts_lib : mts_libs) {
            SelectorNodeOP n (new SelectorNode(mts_lib, i));
            nodes_[0]->add_connection(n);
            nodes_.push_back(n);
        }
    }
    _setup_clash_lists();
}

void
MotifTreeStateSelector::_setup_clash_lists() {
    StringIntMap seen;
    String name, line, clist_name;
    for (auto const & n : nodes_) {
        for (auto const & n2 : n->connections()) {
            if(n->mts_lib()->mtype()  == UNKNOWN) { continue; }
            if(n2->mts_lib()->mtype() == UNKNOWN) { continue; }
            name = resources_path() + "/precomputed/motif_tree_states/";
            name += type_to_str(n->mts_lib()->mtype()) + "_";
            name += type_to_str(n2->mts_lib()->mtype()) + ".clist";
            if(seen.find(name) != seen.end()) { continue; }
            seen[name] = 1;
            if(!file_exists(name)) {
                std::cout << "warning no clash list loaded for " << type_to_str(n->mts_lib()->mtype()) << " and ";
                std::cout << type_to_str(n2->mts_lib()->mtype()) << std::endl;
                continue;
            }
            
            StringIntMapOP clist (new StringIntMap() );
            std::ifstream in;
            in.open(name);
            while ( in.good() ) {
                getline(in, line);
                if(line.size() > 10) { continue; }
                clist->insert( std::pair<String, int>(line, 1));
            }
            in.close();
            clist_name = type_to_str(n->mts_lib()->mtype()) + "-" + type_to_str(n2->mts_lib()->mtype());
            clash_lists_[clist_name] = clist;
            
        }
    }
}


MTSTypePairs const &
MotifTreeStateSelector::get_children_mts(
    MotifTreeStateSearchNodeOP const & node) {
    SelectorNodeOPs connections;
    int lib_type = node->lib_type();
    if(lib_type != -1)         { connections = nodes_[node->lib_type()]->connections(); }
    else                       { connections.push_back(nodes_[0]);  lib_type = 0; }
    
    String clist_name, key;
    StringIntMapOP clist;
    mts_type_pairs_.resize(0);
    for( auto const & c : connections) {
        //if( node->lib_type() != -1 && c->max_uses() <= node->lib_type_usage(c->index()) ) { continue; }
        clist_name = type_to_str(nodes_[lib_type]->mts_lib()->mtype()) + "-" +
                     type_to_str(c->mts_lib()->mtype());
        clist = NULL;
        if(clash_lists_.find(clist_name) != clash_lists_.end() ) {
            clist = clash_lists_[clist_name];
        }
        for (auto const & mts : c->mts_lib()->motif_tree_states()) {
            //key = node->mts()->name() + " " + mts->name();
            //if( clist != NULL && clist->find(key) != clist->end()) { continue; }
            mts_type_pairs_.push_back(MTSTypePair(mts, c->index()));
        }
    }
    return mts_type_pairs_;
}


MotifTreeStateSelector
default_selector(
    MotifTypes types,
    String mode) {
    
    if(types.size() == 0 ) { types.push_back(TWOWAY); }
    
    MotifTreeStateLibraryOPs mts_libs;
    for(auto const & t : types) {
        MotifTreeStateLibraryOP mts_lib ( new MotifTreeStateLibrary(t) );
        mts_libs.push_back(mts_lib);
    }
    
    MotifTreeStateSelector selector(mts_libs, mode);
    return selector;
    
    
}