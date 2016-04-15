//
//  thermo_fluc_simulation_devel.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 11/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//
//

#include "thermo_fluctuation/thermo_fluc_simulation_devel.h"

void
ThermoFlucSimulationDevel::setup_options() {
    options_ = Options();
    options_.add_option(Option("temperature", 298.15f));
    options_.add_option(Option("steps", 100000));
    options_.add_option(Option("record", 0));
    options_.add_option(Option("record_file", "test.out"));
    options_.add_option(Option("cutoff", 5.0f));
    options_.add_option(Option("steric_radius", 2.2f));
    options_.add_option(Option("d_weight", 1.0f));
    options_.add_option(Option("r_weight", 1.0f));
    options_.add_option(Option("record_state", 0));
    options_.add_option(Option("record_all", 0));
    options_.add_option(Option("record_all_file", "test_all.out"));
    options_.add_option(Option("bound_pdb", 0));
    
    update_var_options();
}

void
ThermoFlucSimulationDevel::update_var_options() {
    temperature_    = options_.option<float>("temperature");
    steps_          = options_.option<int>("steps");
    record_         = options_.option<int>("record");
    record_state_   = options_.option<int>("record_state");
    record_all_     = options_.option<int>("record_all");
    record_all_file_= options_.option<String>("record_all_file");
    record_file_    = options_.option<String>("record_file");
    cutoff_         = options_.option<float>("cutoff");
    steric_radius_  = options_.option<float>("steric_radius");
    bound_pdb_      = options_.option<int>("bound_pdb");
    
    std::dynamic_pointer_cast<FrameScorerDevel>(scorer_)->weight_d(options_.option<float>("d_weight"));
    std::dynamic_pointer_cast<FrameScorerDevel>(scorer_)->weight_r(options_.option<float>("r_weight"));

}

void
ThermoFlucSimulationDevel::setup(
    MotifStateEnsembleTreeOP const & mset,
    int ni1,
    int ni2,
    int ei1,
    int ei2) {
    
    sampler_.setup(mset);
    sampler_.temperature(option<float>("temperature"));
    
    ni1_ = ni1;
    ni2_ = ni2;
    ei1_ = ei1;
    ei2_ = ei2;
}


String
ThermoFlucSimulationDevel::static_run() {
    

    float d_weight = options_.option<float>("d_weight");
    float r_weight = options_.option<float>("r_weight");
    
    String results;
    end_state_1_ = sampler_.mst()->get_node(ni1_)->data()->cur_state->end_states()[ei1_];
    end_state_2_ = sampler_.mst()->get_node(ni2_)->data()->cur_state->end_states()[ei2_];
    
    float frame_score_ = end_state_1_->d().distance(end_state_2_->d())*d_weight;
    float r_diff_ = end_state_1_->r().difference(end_state_2_->r());
    end_state_2_->flip();
    float r_diff_flip_ = end_state_1_->r().difference(end_state_2_->r());;
    end_state_2_->flip();
    if(r_diff_ < r_diff_flip_) { r_diff_ = r_diff_flip_; }
    r_diff_ *= r_weight;

    results += std::to_string(frame_score_) + " " + std::to_string(r_diff_) + " " + std::to_string(frame_score_ + r_diff_);
    
    return results;
}

int
ThermoFlucSimulationDevel::run() {
    
    int steps = 0;
    int r = 0;
    int count = 0;
    int clash = 0;
    
    std::cout << record_all_file_ << std::endl;
    
    if(check_nodes_.size() == 0 && check_nodes_2_.size() == 0) {
        check_nodes_ = { 22, 21 };
        check_nodes_2_ = { 1 };
    }
    
    //fixes if length of tecto changes, need to come up with better system!
    
    if(check_nodes_.size() == 2) {
        check_nodes_[0] = sampler_.mst()->last_node()->index();
        check_nodes_[1] = sampler_.mst()->last_node()->index()-1;
    }
        
    std::ofstream out, out_state, out_all;
    if(record_) {
        out.open(record_file_);
        out << "d1,r1,d2,r2,cutoff";
        out << std::endl;
    }
    
    if(record_state_) {
        out_state.open("test_state.out");
        int last = 0;
        int c = 1;
        for(int a = 2; a < sampler_.mst()->size(); a++) {
            if(sampler_.mst()->get_node(a)->data()->cur_state->size() > 4) {
                last = a;
                break;
            }
            out_state << "f" << c << ",";
            c += 1;
        }
        
        c = 1;
        for(int a = last+1; a < sampler_.mst()->size(); a++) {
            out_state << "c" << c << ",";
            c += 1;
        }
        
        out_state << "cutoff" << std::endl;
    }
    
    if(record_all_) {
        out_all.open(record_all_file_);
        int last = 0;
        int c = 1;
        for(int a = 2; a < sampler_.mst()->size(); a++) {
            if(sampler_.mst()->get_node(a)->data()->cur_state->size() > 4) {
                last = a;
                break;
            }
            out_all << "f" << c << "_d," << "f" << c << "_r,";
            c += 1;
        }
        
        c = 1;
        for(int a = last+1; a < sampler_.mst()->size(); a++) {
            out_all << "c" << c << "_d," << "c" << c << "_r,";
            c += 1;
        }
        
        out_all << "cutoff" << std::endl;
    }

    int ncount = 0;
    while (steps < steps_) {
        r = sampler_.next();
        //if(r == 0) { continue; }
        
        clash = 0;
        for(auto const & i : check_nodes_) {
            for(auto const & j : check_nodes_2_) {
                for(auto const & b2 : sampler_.mst()->get_node(i)->data()->cur_state->beads()) {
                    for(auto const & b1 : sampler_.mst()->get_node(j)->data()->cur_state->beads()) {
                        if(b1.distance(b2) < steric_radius_) { clash = 1; }
                    }
                    
                    if(clash) { break; }
                }
                if(clash) { break; }
            }
            if(clash) { break; }
        }
        
        if(clash) {
            //steps++;
            continue;
        }
        
        end_state_1_ = sampler_.mst()->get_node(ni1_)->data()->cur_state->end_states()[ei1_];
        end_state_2_ = sampler_.mst()->get_node(ni2_)->data()->cur_state->end_states()[ei2_];
        
        score_ = scorer_->score(end_state_1_, end_state_2_);
        if(score_ < cutoff_) {            
            
            count += 1;
            
            if(bound_pdb_) {
                sampler_.mst()->to_motif_tree()->to_pdb("bound.pdb");
                exit(0);
            }
        }
        
        if(record_) {
            out << vector_to_str(end_state_1_->d()) << "," << matrix_to_str(end_state_1_->r()) << "," <<vector_to_str(end_state_2_->d()) << "," << matrix_to_str(end_state_2_->r()) << "," << cutoff_;
            out << std::endl;
        }
        
        if(record_state_) {
            int last = 0;
            for(int a = 2; a < sampler_.mst()->size(); a++) {
                if(sampler_.mst()->get_node(a)->data()->cur_state->size() > 4) {
                    last = a;
                    break;
                }
                out_state << sampler_.mst()->get_node(a)->data()->cur_state->name() <<  ",";
            }
            
            for(int a = last+1; a < sampler_.mst()->size(); a++) {
                if(sampler_.mst()->get_node(a)->data()->cur_state->size() > 4) {
                    last = a;
                    break;
                }
                out_state << sampler_.mst()->get_node(a)->data()->cur_state->name() <<  ",";
            }

            out_state << score_ << std::endl;
        }
        
        if(record_all_) {
            int last = 0;
            for(int a = 2; a < sampler_.mst()->size(); a++) {
                if(sampler_.mst()->get_node(a)->data()->cur_state->size() > 4) {
                    last = a;
                    break;
                }
                out_all << vector_to_str(sampler_.mst()->get_node(a)->data()->cur_state->end_states()[1]->d()) <<  ",";
                out_all << matrix_to_str(sampler_.mst()->get_node(a)->data()->cur_state->end_states()[1]->r()) <<  ",";
            }
            
            for(int a = last+1; a < sampler_.mst()->size(); a++) {
                if(sampler_.mst()->get_node(a)->data()->cur_state->size() > 4) {
                    last = a;
                    break;
                }
                out_all << vector_to_str(sampler_.mst()->get_node(a)->data()->cur_state->end_states()[1]->d()) <<  ",";
                out_all << matrix_to_str(sampler_.mst()->get_node(a)->data()->cur_state->end_states()[1]->r()) <<  ",";            }
            
            out_all << score_ << std::endl;
        }
        
        steps++;
    }
    
    if(record_) {
        out.close();
    }
    
    if(record_state_) {
        out_state.close();
    }
    
    if(record_all_) {
        out_all.close();
    }
    
    return count;
}