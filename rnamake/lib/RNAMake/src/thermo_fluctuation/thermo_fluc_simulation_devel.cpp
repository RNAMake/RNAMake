//
//  thermo_fluc_simulation_devel.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 11/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//
//

#include "thermo_fluctuation/thermo_fluc_simulation_devel.h"
#include <sys/stat.h>


void
ThermoFlucSimulationDevel::setup(
    MotifStateEnsembleTreeOP const & mset,
    int ni1,
    int ni2,
    int ei1,
    int ei2) {
    
    sampler_.setup(mset);
    sampler_.temperature(get_float_option("temperature"));
    
    ni1_ = ni1;
    ni2_ = ni2;
    ei1_ = ei1;
    ei2_ = ei2;
    
    // make sure node information that is given is correct
    try {
        sampler_.mst()->get_node(ni1_);
    }
    catch(TreeException const & e) {
        throw ThermoFlucSimulationException(
            "cannot setup thermo fluc simulation as there is no node with index: " +
            std::to_string(ni1_) + " in given tree");
    }
    
    try {
        sampler_.mst()->get_node(ni2_);
    }
    catch(TreeException const & e) {
        throw ThermoFlucSimulationException(
            "cannot setup thermo fluc simulation as there is no node with index: " +
            std::to_string(ni2_) + " in given tree");
    }
    
    if(sampler_.mst()->get_node(ni1_)->data()->cur_state->end_states().size() <= ei1_) {
        throw ThermoFlucSimulationException(
            "cannot setup thermo fluc simulation as node: " + std::to_string(ni1_) +
            " does not have a " + std::to_string(ei1_) + " end");
    }
    
    if(sampler_.mst()->get_node(ni2_)->data()->cur_state->end_states().size() <= ei2_) {
        throw ThermoFlucSimulationException(
            "cannot setup thermo fluc simulation as node: " + std::to_string(ni2_) +
            " does not have a " + std::to_string(ei2_) + " end");
    }
    
    //parse sterics information
    if(get_string_option("steric_nodes") == "") { return; }
    auto steric_node_str = get_string_option("steric_nodes");
    auto spl = split_str_by_delimiter(steric_node_str, ":");
    if(spl.size() != 2) {
        throw ThermoFlucSimulationException(
            "incorrect format for steric_nodes option, must be in the form NodeSet1:NodeSet2,"
            " example: '22,21:1' which check sterics of nodes 22 and 21 against node 1");
    }
    
    auto node_set_1 = split_str_by_delimiter(spl[0], ",");
    for(auto const & n_num : node_set_1) {
        auto n_index = std::stoi(n_num);
        
        try {
            sampler_.mst()->get_node(n_index);
        }
        catch(TreeException const & e) {
            throw ThermoFlucSimulationException(
                "cannot setup sterics for simulation, specified node: " + std::to_string(n_index) +
                " does not exist");
        }
        
        check_nodes_1_.push_back(n_index);
    }
    
    auto node_set_2 = split_str_by_delimiter(spl[1], ",");
    for(auto const & n_num : node_set_2) {
        auto n_index = std::stoi(n_num);
        
        try {
            sampler_.mst()->get_node(n_index);
        }
        catch(TreeException const & e) {
            throw ThermoFlucSimulationException(
                "cannot setup sterics for simulation, specified node: " + std::to_string(n_index) +
                " does not exist");
        }
        
        check_nodes_2_.push_back(n_index);
    }
    

}


String
ThermoFlucSimulationDevel::static_run() {
    

    /*float d_weight = options_.option<float>("d_weight");
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

    results += std::to_string(frame_score_) + " " + std::to_string(r_diff_) + " " + std::to_string(frame_score_ + r_diff_);*/
    
    //return results;
    return String("");
}


int
ThermoFlucSimulationDevel::run() {
    
    int steps = 0;
    int r = 0;
    int count = 0;
    int clash = 0;
    int pdb_count = 0;

    auto r1_trans = Matrix();
    auto r2_trans = Matrix();
    auto r_result = Matrix();
    auto rot      = Matrix();

    auto mf = MotifFactory();

    std::ofstream out, out_state, out_all, out_motifs, out_dump_state;

    if(dump_state_) {
        out_dump_state.open("state.out");
    }

    if(record_) {
        // default setup
        if(logger_ == nullptr) { logger_ = std::make_shared<ThermoFlucSimulationLogger>(); }

        logger_->setup(sampler_.mst(), ni1_, ni2_, ei1_, ei2_);
    }

    while (steps < steps_) {
        r = sampler_.next();
        //TODO revisit if moving to the next state gives better results
        if(r == 0) { continue; }

        clash = _check_sterics();
        //clash = 0;
        
        if(clash) {
            steps++;
            continue;
        }

        end_state_1_ = sampler_.mst()->get_node(ni1_)->data()->get_end_state(ei1_);
        end_state_2_ = sampler_.mst()->get_node(ni2_)->data()->get_end_state(ei2_);
        
        score_ = scorer_->score(end_state_1_, end_state_2_);
        if(score_ < cutoff_) { count += 1; }

        if(record_) {
            if     (record_only_bound_ && score_ <= cutoff_)      { logger_->log(sampler_.mst(), score_); }
            else if(record_only_unbound_ && score_ > cutoff_)     { logger_->log(sampler_.mst(), score_); }
            else if(!record_only_bound_ && !record_only_unbound_) { logger_->log(sampler_.mst(), score_); }
        }
        if(dump_state_) {
            if     (record_only_bound_ && score_ <= cutoff_)      {
                out_dump_state << sampler_.mst()->to_motif_tree()->to_str() << std::endl;
            }
            else if(record_only_unbound_ && score_ > cutoff_)     {
                out_dump_state << sampler_.mst()->to_motif_tree()->to_str() << std::endl;
            }
            else if(!record_only_bound_ && !record_only_unbound_) {
                out_dump_state << sampler_.mst()->to_motif_tree()->to_str() << std::endl;
            }
        }
        if(dump_pdbs_) {
            if     (record_only_bound_ && score_ <= cutoff_)      {
                sampler_.mst()->to_motif_tree()->to_pdb("test."+std::to_string(pdb_count)+".pdb", 1, 1);
                pdb_count += 1;
            }
            else if(record_only_unbound_ && score_ > cutoff_)     {
                sampler_.mst()->to_motif_tree()->to_pdb("test."+std::to_string(pdb_count)+".pdb", 1, 1);
                pdb_count += 1;
            }
            else if(!record_only_bound_ && !record_only_unbound_) {
                sampler_.mst()->to_motif_tree()->to_pdb("test."+std::to_string(pdb_count)+".pdb", 1, 1);
                pdb_count += 1;
            }

        }

        steps++;
    }
    if(record_) {
        logger_->finalize();
    }
    return count;
}


void
ThermoFlucSimulationDevel::setup_options() {
    options_.add_option("temperature", 298.15f, OptionType::FLOAT);
    options_.add_option("steps", 100000, OptionType::INT);
    options_.add_option("cutoff", 4.5f, OptionType::FLOAT);
    options_.add_option("steric_nodes", "", OptionType::STRING);
    options_.add_option("record", false, OptionType::BOOL);
    options_.add_option("record_file", "test.out", OptionType::STRING);
    options_.add_option("record_state", false, OptionType::BOOL);
    options_.add_option("record_all", false, OptionType::BOOL);
    options_.add_option("record_all_file", "test_all.out", OptionType::STRING);
    options_.add_option("unbound_pdbs", false, OptionType::BOOL);
    options_.add_option("dump_state", false, OptionType::BOOL);
    options_.add_option("dump_pdbs", false, OptionType::BOOL);
    options_.add_option("record_only_bound", false, OptionType::BOOL);
    options_.add_option("record_only_unbound", false, OptionType::BOOL);
    options_.add_option("steric_radius", 2.2f, OptionType::FLOAT);
    options_.lock_option_adding();
    
    /*
     options_.add_option(Option("d_weight", 1.0f));
     options_.add_option(Option("r_weight", 1.0f));
     options_.add_option(Option("record_state", 0));
     options_.add_option(Option("record_all", 0));
     options_.add_option(Option("record_all_file", "test_all.out"));
     options_.add_option(Option("bound_pdb", 0));*/
    
    update_var_options();
}

void
ThermoFlucSimulationDevel::update_var_options() {
    temperature_          = options_.get_float("temperature");
    steps_                = options_.get_int("steps");
    cutoff_               = options_.get_float("cutoff");
    record_               = options_.get_bool("record");
    record_all_           = options_.get_bool("record_all");
    record_all_file_      = options_.get_string("record_all_file");
    record_file_          = options_.get_string("record_file");
    record_state_         = options_.get_bool("record_state");
    steric_radius_        = options_.get_float("steric_radius");
    dump_state_           = options_.get_bool("dump_state");
    dump_pdbs_            = options_.get_bool("dump_pdbs");
    record_only_bound_    = options_.get_bool("record_only_bound");
    record_only_unbound_  = options_.get_bool("record_only_unbound");

}
