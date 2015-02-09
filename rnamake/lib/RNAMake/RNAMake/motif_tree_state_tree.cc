//
//  motif_tree_state_tree.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <queue>
#include "motif_tree_state_tree.h"
#include "motif_tree.h"
#include "motif.h"
#include "basepair.h"

MotifTreeStateTree::MotifTreeStateTree() {
    clash_radius_ = 2.9;
    sterics_ = 1;
    nodes_ = MotifTreeStateNodeOPs();
    aligner_ = MotifTreeStateNodeAligner();
    MotifTreeStateNodeOP head ( new MotifTreeStateNode(MotifTreeStateOP(new MotifTreeState(ref_mts())), 0, NULL, 0));
    nodes_.push_back(head);
    last_node_ = head;
}

MotifTreeStateTree::MotifTreeStateTree(MotifTreeStateOP const & mts) {
    clash_radius_ = 2.9;
    sterics_ = 1;
    nodes_ = MotifTreeStateNodeOPs();
    aligner_ = MotifTreeStateNodeAligner();
    MotifTreeStateNodeOP head ( new MotifTreeStateNode(mts, 0, NULL, 0));
    nodes_.push_back(head);
    last_node_ = head;

}

MotifTreeStateNodeOP
MotifTreeStateTree::add_state(
    MotifTreeStateOP const & mts,
    MotifTreeStateNodeOP const & cparent,
    int const cparent_end) {
    
    MotifTreeStateNodeOP parent = cparent;
    Ints indices;
    if(parent == NULL) { parent = last_node_; }
    if(cparent_end == -1) {
        indices = parent->available_ends();
    }
    else {
        indices.push_back(cparent_end);
    }
    
    MotifTreeStateNodeOP new_node ( new MotifTreeStateNode(mts, (int)nodes_.size(), parent, 0));
    int success = 0;

    for (auto const & i : indices) {
        BasepairStateOP state = parent->states()[i];
        aligner_.transform_state(state, parent, new_node);
        aligner_.transform_beads(new_node);
        if(sterics_ == 1) {
            if(_steric_clash(new_node)) {
                continue;
            }
        }
        parent->add_child(new_node, i);
        success=1;
        break;
    }
    
    if(!success) { return NULL; }
    
    nodes_.push_back(new_node);
    last_node_ = new_node;
    return new_node;
    
}

int
MotifTreeStateTree::_steric_clash(
    MotifTreeStateNodeOP const & new_node) {
    float dist;
    for( auto const & n : nodes_) {
        for (auto const & b1 : new_node->beads()) {
            for (auto const & b2 : n->beads()) {
                dist = b1.distance(b2);
                if(dist < clash_radius_) { return 1; }
            }
        }
    }
    return 0;
    
    
}

MotifTree
MotifTreeStateTree::to_motiftree() {
    int i = -1;
    MotifTree mt;
    ResidueTypeSet rts;
    for (auto const & n : nodes_) {
        i++;
        if(i == 0) {
            if(n->mts()->name().compare("start") == 0) {
                mt = MotifTree();
                mt.sterics(0);
            }
            else {
                MotifOP m ( new Motif( n->mts()->build_string(), rts ));
                mt = MotifTree(m);
                //mt.sterics(0);
            }
            continue;
        }
        MotifOP m ( new Motif( n->mts()->build_string(), rts ));
        int pos = (int)(std::find(nodes_.begin(), nodes_.end(), n->parent()) - nodes_.begin());
        MotifTreeNodeOP parent = mt.nodes()[pos];
        int parent_index = n->parent_end_index();
        //TODO: not sure why including flip screws stuff up?
        MotifTreeNodeOP mt_node = mt.add_motif(m, parent, n->mts()->start_index(), parent_index, n->mts()->flip());
        //MotifTreeNodeOP mt_node = mt.add_motif(m, parent, n->mts()->start_index(), parent_index);
        if(mt_node == NULL) {
            std::cout << i << " " << n->mts()->name() << " " << n->mts()->flip() << std::endl;
            mt.sterics(0);
            MotifTreeNodeOP mt_node = mt.add_motif(m, parent, n->mts()->start_index(), parent_index, n->mts()->flip());
            mt.write_pdbs();
            exit(0);
        }
    }
    return mt;
}

int
MotifTreeStateTree::replace_state(
    MotifTreeStateNodeOP const & node,
    MotifTreeStateOP const & mts) {
    
    //check node topology bad things happen if it changes
    if( node->mts()->end_states().size() != mts->end_states().size()) { return 0; }
    for (int i = 0; i < node->mts()->end_states().size(); i++) {
        if(node->mts()->end_states()[i] == NULL && mts->end_states()[i] != NULL ) { return 0; }
        if(node->mts()->end_states()[i] != NULL && mts->end_states()[i] == NULL ) { return 0; }
    }
    
    MotifTreeStateOP old_mts = node->mts();
    node->replace_mts(mts);
    std::queue<MotifTreeStateNodeOP> open_nodes;
    MotifTreeStateNodeOPs seen_nodes;
    MotifTreeStateNodeOP current, parent;
    BasepairStateOP parent_end;
    open_nodes.push(node);
    while (! open_nodes.empty() ) {
        current = open_nodes.front();
        open_nodes.pop();
        seen_nodes.push_back(current);
        parent = current->parent();
        if(parent == NULL) { continue; }
        parent_end = current->parent_end();
        aligner_.transform_state(parent_end, parent, current);
        aligner_.transform_beads(current);
        for (auto const & c : current->children()) {
            open_nodes.push(c);
        }
    }
    
    //To not check sterics
    //return 1;
    
    int clash = 0;
    float dist = 0;
    for(int i = 0; i < nodes_.size(); i++) {
        for(int j = i+1; j < nodes_.size(); j++) {
            for( auto const & b1: nodes_[i]->beads()) {
                for (auto const & b2 : nodes_[j]->beads()) {
                    dist = b1.distance(b2);
                    if(dist < clash_radius_) { clash = 1; }
                }
            }
        }
    }
    
    if(!clash) { return 1; }
    
    node->replace_mts(old_mts);
    
    for (auto const & n : seen_nodes) {
        parent = n->parent();
        parent_end = n->parent_end();
        if(parent == NULL) { continue; }
        aligner_.transform_state(parent_end, parent, current);
        aligner_.transform_beads(current);
    }
    
    return 0;
}



MotifTreeState
ref_mts() {
    Motif m = ref_motif();
    BasepairOP ref_bp = m.ends()[0];
    Points beads;
    for( auto const & b : ref_bp->res1()->get_beads()) {
        if(b.btype() != PHOS ) {
            beads.push_back(b.center());
        }
    }
    for( auto const & b : ref_bp->res2()->get_beads()) {
        if(b.btype() != PHOS ) {
            beads.push_back(b.center());
        }
    }
    BasepairStateOPs states;
    states.push_back( BasepairStateOP(new BasepairState(ref_bp->state())));
    MotifTreeState start_mts ("start", 1, 0, 0, beads, states, 0, ref_motif().to_str());
    return start_mts;
}

