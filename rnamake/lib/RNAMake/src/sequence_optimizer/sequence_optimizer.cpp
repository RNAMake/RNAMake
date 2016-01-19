//
//  sequence_optimizer.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/11/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "motif_data_structures/motif_topology.h"
#include "sequence_optimizer/sequence_optimizer.h"

SequenceOptimizerResultOP
SequenceOptimizer::optimize(
    MotifGraphOP & mg,
    int node_i,
    int node_j,
    int end_i,
    int end_j) {
    
    auto dss = mg->designable_secondary_structure();
    auto org_dss = mg->designable_secondary_structure();
    auto node_name = mg->get_node(node_j)->data()->name();
    int count = 0;
    for (auto const & n : *mg) {
        if(n->data()->name() == node_name) { count++; }
    }
    assert(count == 1 && "cannot optimize sequence too many nodes of the same name");
    
    designer_ = eternabot::SequenceDesigner();
    designer_.setup();
    
    designer_results_ = designer_.design(dss);
    int i = 0;
    int new_node_j = 0;
    Point d1, d2;
    float dist;
    
    auto best_seq = String();
    auto best = 100000;
    
    for(auto const & r : designer_results_) {
        i++;

        if(r->sequence.length() < 2) { break; }
        
        dss->replace_sequence(r->sequence);

        
        mg->replace_helical_sequence(dss);
        mg->write_pdbs();
      
        
        mt_ = graph_to_tree(mg, mg->get_node(node_i),
                            mg->get_node(node_j)->data()->ends()[end_j]);
        
        
        new_node_j = 0;
        for(auto const & n : *mt_) {
            if(n->data()->name() == node_name) {
                new_node_j = n->index();
                break;
            }
        }
        
        d1 = mt_->get_node(new_node_j)->data()->ends()[end_j]->d();
        d2 = mt_->last_node()->data()->ends()[1]->d();
        dist = d1.distance(d2);

        if(best > dist) {
            best = dist;
            best_seq = r->sequence;
        }
        
        if(i > 10) {
            break;
        }
        
    }
    
    dss->replace_sequence(best_seq);
    mg->replace_helical_sequence(dss);
    mg->write_pdbs();
    exit(0);
    mt_ = graph_to_tree(mg, mg->get_node(node_i),
                        mg->get_node(node_j)->data()->ends()[end_j]);
    
    return std::make_shared<SequenceOptimizerResult>(mt_, best);

    
    
}












