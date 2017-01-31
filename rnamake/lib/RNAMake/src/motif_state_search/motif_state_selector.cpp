//
//  motif_state_selector.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_state_selector.h"


void
MotifStateSelector::connect(
    String const & name_i,
    String const & name_j) {
    
    int i = -1, j = -1;
    for(auto const & n : graph_) {
        if(n->data()->name == name_i) { i = n->index(); }
        if(n->data()->name == name_j) { j = n->index(); }
    }
    
    if(i == -1 || j == -1) {
        throw std::runtime_error("could not connect nodes in MotifStateSelector::connect");
    }
    graph_.connect(i, j);
}

MotifStateandTypes const &
MotifStateSelector::get_children_ms(
    MotifStateSearchNodeOP const & node) {

    motif_states_and_types_.reset();
    //beginning of search start on first node
    if(node->ntype() == -1) {
        motif_states_and_types_.need_resize((int)graph_.get_node(0)->data()->motif_states.size());
        for(auto const & ms : graph_.get_node(0)->data()->motif_states) {
            motif_states_and_types_.add(ms, 0);
        }
    }
    
    else {
        for(auto const & c : graph_.get_node(node->ntype())->connections()) {
            auto partner =  c->partner(graph_.get_node(node->ntype())->index());
            //did use all that we needed of this type
            if(partner->data()->max_uses <= node->node_type_usage(partner->index())) {
                continue;
            }
            motif_states_and_types_.need_resize((int)partner->data()->motif_states.size());
            for(auto const & ms : partner->data()->motif_states) {
                motif_states_and_types_.add(ms, partner->index());
            }
        }
    }
    
    return motif_states_and_types_;
}

int
MotifStateSelector::is_valid_solution(
    MotifStateSearchNodeOP const & current) {
    int i = 0;
    for(auto const & n : graph_) {
        if(n->data()->required_uses > current->node_type_usage(i)) {
            return 0;
        }
        i++;
    }
    return 1;
}

float
MotifStateSelector::score(
    MotifStateSearchNodeOP const & current) {
    float diff = 0;
    int i = 0;
    for(auto const & n : graph_) {
        if(n->data()->required_uses > current->node_type_usage(i)) {
            diff += (n->data()->required_uses - current->node_type_usage(i));
        }
    }
    
    return diff;
}

MotifStateSelectorOP
default_selector() {
    auto selector = std::make_shared<MotifStateSelector>(MSS_HelixFlank());
    //selector->add("unique_twoway");
    selector->add("twoway");
    
    
    selector = std::make_shared<MotifStateSelector>(MotifStateSelector());
    selector->add("twoway");
    selector->add("ideal_helices_min");
    selector->connect("twoway", "ideal_helices_min");
    
    return selector;
}
