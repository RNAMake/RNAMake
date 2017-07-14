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

void
MotifStateSelector::connect(
        int i,
        int j) {
    graph_.connect(i, j);
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
    selector->add("unique_twoway");
    //selector->add("twoway");


    //selector = std::make_shared<MotifStateSelector>(MotifStateSelector());
    //selector->add("twoway");
    //selector->add("ideal_helices_min");selector->add("nway");
    //selector->connect("twoway", "ideal_helices_min");

    return selector;
}


