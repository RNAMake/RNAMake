//
//  sequence_optimization.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/11/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure/secondary_structure_parser.h"
#include "motif_data_structure/motif_topology.h"
#include "sequence_optimization/sequence_optimizer.h"

namespace sequence_optimization {

/*
void
SequenceOptimizer::setup_options() {
    options_.add_option("sub_sequence", String(""), base::OptionType::STRING);
    options_.add_option("solutions", 100, base::OptionType::INT);
    options_.lock_option_adding();
    update_var_options();
}

void
SequenceOptimizer::update_var_options() {
    solutions_ = get_int_option("solutions");
    
}


OptimizedSequenceOPs
SequenceOptimizer::get_optimized_sequences(
    motif_data_structure::MotifGraphOP & mg,
    util::Uuid const & uuid_1,
    util::Uuid const & uuid_2,
    int end_1,
    int end_2) {
    
    auto unaligned_nodes = mg->unaligned_nodes();
    if(unaligned_nodes.size() != 1) {
        throw std::runtime_error("unsure where to base motif tree");
    }
        
    auto sols = OptimizedSequenceOPs();
    auto dss = mg->designable_secondary_structure();
    auto start_dss = mg->designable_secondary_structure();

    
    designer_ = eternabot::SequenceDesigner();
    designer_.setup();
    
    //hack to only optimize a part of the sequence at a time, so vienna doesn't get confused
    //with multiple strands
    if(get_string_option("sub_sequence") != "") {
        sub_sequence_ = true;
        auto spl = base::split_str_by_delimiter(get_string_option("sub_sequence"), ":");
        start_ = std::stoi(spl[0]);
        end_   = std::stoi(spl[1]);
        if (end_ == -1) {
            end_ = start_dss->sequence().length();
        }
        
        auto seq = start_dss->sequence().substr(start_, (end_-start_));
        auto db  = start_dss->dot_bracket().substr(start_, (end_-start_));
        auto parser = secondary_structure::Parser();
        auto new_p = parser.parse_to_pose(seq, db);
        start_dss = new_p;
    }
    
    designer_results_ = designer_.design(start_dss);
    auto org_seq = dss->sequence();
    auto designed_seq = String("");
    int i = 0;
    int new_node_j = 0;
    math::Point d1, d2;
    float dist;
    
    auto best_seq = String();
    auto best = 100000;
    
    auto end_state_1 = structure::BasepairStateOP();
    auto end_state_2 = structure::BasepairStateOP();
    auto end_state_2_flip = structure::BasepairStateOP();
    
    auto node_1 = mg->get_node(uuid_1);
    
    for(auto const & r : designer_results_) {
        i++;
        
    
        if(r->sequence.length() < 2) { break; }
        designed_seq = get_final_sequence(r->sequence, org_seq);
        
        dss->replace_sequence(designed_seq);
        mg->replace_helical_sequence(dss);
        
        mt_ = graph_to_tree(mg, unaligned_nodes[0],
                            node_1->data()->ends()[end_1]);
        
        
        end_state_1 = mt_->get_node(uuid_1)->data()->ends()[end_1]->state();
        end_state_2 = mt_->get_node(uuid_2)->data()->ends()[end_2]->state();
        end_state_2->flip();
        end_state_2_flip = std::make_shared<structure::BasepairState>(end_state_2->copy());
        end_state_2->flip();
        dist = new_score_function_new(end_state_1, end_state_2, end_state_2_flip);
        
        sols.push_back(std::make_shared<OptimizedSequence>(designed_seq, dist, r->score));
        if(i > solutions_-1) {
            break;
        }
        
    }
    
    return sols;
}


String
SequenceOptimizer::get_final_sequence(
    String const & designed_seq,
    String const & org_seq) {
    
    if(!sub_sequence_) { return designed_seq; }
    
    String full_seq = "";
    if(start_ > 0) {
        full_seq += org_seq.substr(0, start_);
    }
    full_seq += designed_seq;
    if(end_ < org_seq.length()) {
        full_seq += org_seq.substr(end_);
    }
    return full_seq;
    
}
*/

}






