//
//  sequence_optimizer.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/11/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure/secondary_structure_parser.h"
#include "motif_data_structures/motif_topology.h"
#include "sequence_optimizer/sequence_optimizer.h"

void
SequenceOptimizer::setup_options() {
    options_.add_option("sub_sequence", String(""), OptionType::STRING);
    options_.add_option("solutions", 100, OptionType::INT);
    options_.lock_option_adding();
    update_var_options();
}

void
SequenceOptimizer::update_var_options() {
    solutions_ = get_int_option("solutions");
    
}


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

        std::cout << dist << " " << r->score << " " << r->sequence << std::endl;
   
        
        if(best > dist) {
            best = dist;
            best_seq = r->sequence;
        }
        
        if(i > 100) {
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


OptimizedSequenceOPs
SequenceOptimizer::get_optimized_sequences(
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
    
    auto sols = OptimizedSequenceOPs();
    auto end_state_1 = BasepairStateOP();
    auto end_state_2 = BasepairStateOP();
    auto end_state_2_flip = BasepairStateOP();
    
    for(auto const & r : designer_results_) {
        i++;
        
        if(r->sequence.length() < 2) { break; }
        
        
        dss->replace_sequence(r->sequence);
 
        mt_ = graph_to_tree(mg, mg->get_node(node_i),
                            mg->get_node(node_j)->data()->ends()[end_j]);
        
        new_node_j = 0;
        for(auto const & n : *mt_) {
            if(n->data()->name() == node_name) {
                new_node_j = n->index();
                break;
            }
        }
        
        end_state_1 = mt_->get_node(new_node_j)->data()->ends()[end_j]->state();
        exit(0);
        
        end_state_2 = mt_->get_node(new_node_j)->children()[end_j]->data()->ends()[1]->state();
        end_state_2->flip();
        end_state_2_flip = std::make_shared<BasepairState>(end_state_2->copy());
        end_state_2->flip();
        dist = new_score_function_new(end_state_1, end_state_2, end_state_2_flip);
        
        std::cout << " " << i << " " << dist << std::endl;
        mg->to_pdb("test."+std::to_string(i) + ".pdb");

        sols.push_back(std::make_shared<OptimizedSequence>(r->sequence, dist, r->score));
     
        if(i > 10) {
            break;
        }
        
    }
    
    return sols;
}


OptimizedSequenceOPs
SequenceOptimizer::get_optimized_sequences_2(
    MotifGraphOP & mg,
    Uuid const & uuid_1,
    Uuid const & uuid_2,
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
        auto spl = split_str_by_delimiter(get_string_option("sub_sequence"), ":");
        start_ = std::stoi(spl[0]);
        end_   = std::stoi(spl[1]);
        if (end_ == -1) {
            end_ = start_dss->sequence().length();
        }
        
        auto seq = start_dss->sequence().substr(start_, (end_-start_));
        auto db  = start_dss->dot_bracket().substr(start_, (end_-start_));
        auto parser = sstruct::SecondaryStructureParser();
        auto new_p = parser.parse_to_pose(seq, db);
        start_dss = new_p;
    }
    
    designer_results_ = designer_.design(start_dss);
    auto org_seq = dss->sequence();
    auto designed_seq = String("");
    int i = 0;
    int new_node_j = 0;
    Point d1, d2;
    float dist;
    
    auto best_seq = String();
    auto best = 100000;
    
    auto end_state_1 = BasepairStateOP();
    auto end_state_2 = BasepairStateOP();
    auto end_state_2_flip = BasepairStateOP();
    
    auto node_1 = mg->get_node_by_id(uuid_1);
    
    for(auto const & r : designer_results_) {
        i++;
        
    
        if(r->sequence.length() < 2) { break; }
        designed_seq = get_final_sequence(r->sequence, org_seq);
        
        dss->replace_sequence(designed_seq);
        mg->replace_helical_sequence(dss);
        
        mt_ = graph_to_tree(mg, unaligned_nodes[0],
                            node_1->data()->ends()[end_1]);
        
        
        end_state_1 = mt_->get_node_by_id(uuid_1)->data()->ends()[end_1]->state();
        end_state_2 = mt_->get_node_by_id(uuid_2)->data()->ends()[end_2]->state();
        end_state_2->flip();
        end_state_2_flip = std::make_shared<BasepairState>(end_state_2->copy());
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









