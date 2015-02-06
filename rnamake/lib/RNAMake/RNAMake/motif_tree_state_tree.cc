//
//  motif_tree_state_tree.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_state_tree.h"
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

MotifTreeStateNodeOP
MotifTreeStateTree::add_state(
    MotifTreeStateOP const & mts,
    MotifTreeStateNodeOP const & cparent,
    BasepairStateOP const & parent_end) {
    
    MotifTreeStateNodeOP parent = cparent;
    Ints indices;
    if(parent == NULL) { parent = last_node_; }
    indices = parent->available_ends();
    
    MotifTreeStateNodeOP new_node ( new MotifTreeStateNode(mts, (int)nodes_.size(), parent, 0));
    int success = 0;
    for (auto const & i : indices) {
        BasepairStateOP state = parent->states()[i];
        aligner_.transform_state(state, parent, new_node);
        //aligner_.transform_beads(new_node);
        if(sterics_ == 1) {
        //    if(_steric_clash(new_node)) {
        //        continue;
        //    }
        }
        //std::cout << parent->states().size() << std::endl;
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

