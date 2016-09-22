//
//  thermo_fluc_simulation.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "thermo_fluctuation/thermo_fluc_simulation.h"

ThermoFlucSimulation::ThermoFlucSimulation() {
    scorer_ = std::make_shared<FrameScorer>(FrameScorer());
    sampler_ = ThermoFlucSampler();
    options_ = Options();
    check_nodes_1_ = Ints{};
    check_nodes_2_ = Ints{};
    setup_ = 0;
    setup_options();

}


void
ThermoFlucSimulation::setup(
    MotifStateEnsembleTreeOP const & mset,
    int ni1,
    int ni2,
    int ei1,
    int ei2) {
    
    setup_ = 1;
    sampler_.setup(mset);
    sampler_.temperature(temperature_);
    
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

int
ThermoFlucSimulation::run() {
    
    if(!setup_) {
        throw ThermoFlucSimulationException(
            "tried to execute run() without first calling setup, this will segfault");
    }
    
    int steps = 0;
    int r = 0;
    int count = 0;
    int clash = 0;
    //Ints check_nodes = { 22, 21 };
    //Ints check_nodes_2 = { 1 };
    
    while (steps < steps_) {
        //if(r == 0) { continue; }
        
        r = sampler_.next();

        clash = _check_sterics();
        if(clash) { steps++; continue; }
        
        end_state_1_ = sampler_.mst()->get_node(ni1_)->data()->get_end_state(ei1_);
        end_state_2_ = sampler_.mst()->get_node(ni2_)->data()->get_end_state(ei2_);

        score_ = scorer_->score(end_state_1_, end_state_2_);
        
        if(score_ < cutoff_) { count += 1; }
        
        steps++;
    }
    
    return count;
}


//option functions /////////////////////////////////////////////////////////////////////////////////


void
ThermoFlucSimulation::setup_options() {
    options_.add_option("temperature", 298.15f, OptionType::FLOAT);
    options_.add_option("steps", 100000, OptionType::INT);
    options_.add_option("record", false, OptionType::BOOL);
    options_.add_option("cutoff", 4.5f, OptionType::FLOAT);
    options_.add_option("steric_nodes", "", OptionType::STRING);
    options_.lock_option_adding();
    update_var_options();
}

void
ThermoFlucSimulation::update_var_options() {
    temperature_   = options_.get_float("temperature");
    steps_         = options_.get_int("steps");
    cutoff_        = options_.get_float("cutoff");
    record_        = options_.get_bool("record");
}
