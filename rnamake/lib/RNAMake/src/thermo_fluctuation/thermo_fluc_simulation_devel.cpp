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

    std::ofstream out, out_state, out_all, out_motifs;

    if(bound_pdbs_) {
        out_motifs.open("motifs.out");
    }
    if(record_) {
        out.open(record_file_);
        out << "d1,r1,d2,r2,cutoff,score";
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
        out_all << "fstart_d,fstart_r,";
        for(int a = 2; a < sampler_.mst()->size(); a++) {
            if(sampler_.mst()->get_node(a)->data()->cur_state->size() > 4) {
                last = a;
                break;
            }
            out_all << "f" << c << "_d," << "f" << c << "_r,";
            c += 1;
        }
        
        c = 1;
        out_all << "cstart_d,cstart_r,";
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

        clash = _check_sterics();
        //clash = 0;
        
        if(clash) {
            steps++;
            continue;
        }
        
        end_state_1_ = sampler_.mst()->get_node(ni1_)->data()->get_end_state(ei1_);
        end_state_2_ = sampler_.mst()->get_node(ni2_)->data()->get_end_state(ei2_);
        
        score_ = scorer_->score(end_state_1_, end_state_2_);
        if(score_ < cutoff_) {
            
            count += 1;

            if(bound_pdbs_ && !all_pdbs_) {
                try {
                    //sampler_.mst()->to_motif_tree()->to_pdb("bound." + std::to_string(pdb_count) + ".pdb", 1);
                    auto mt = sampler_.mst()->to_motif_tree();
                    auto bp1 = mt->get_node(ni1_)->data()->ends()[ei1_];
                    auto bp2 = mt->get_node(ni2_)->data()->ends()[ei2_];
                    auto m = mf.motif_from_bps(BasepairOPs{bp1, bp2});
                    //m->to_pdb("bound." + std::to_string(pdb_count) + ".pdb", 1);
                    out_motifs << m->to_str() << std::endl;
                    pdb_count += 1;
                }
                catch(...) {}
            }

            /*if(bound_pdb_ > bound_pdb_count) {
                try {
                    const int dir_err = mkdir(String("nodes_" + std::to_string(bound_pdb_count)).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                    sampler_.mst()->write_pdbs("nodes_" + std::to_string(bound_pdb_count) + "/nodes");
                 
                    //sampler_.mst()->to_motif_tree()->to_pdb("bound."+std::to_string(bound_pdb_count)+
                     //                                       ".pdb");
                    bound_pdb_count++;
                }
                catch(...) {}
            }*/
        }

        else {
            if(unbound_pdbs_ && !all_pdbs_) {
                try {
                    sampler_.mst()->to_motif_tree()->to_pdb("unbound." + std::to_string(pdb_count) + ".pdb", 1);
                    pdb_count += 1;
                }
                catch(...) {}
            }
        }

        if(all_pdbs_) {
            sampler_.mst()->to_motif_tree()->to_pdb("all." + std::to_string(pdb_count) + ".pdb", 1);
            pdb_count += 1;
        }
        
        if(record_) {

            out << vector_to_str(end_state_1_->d()) << "," << matrix_to_str(end_state_1_->r()) << "," <<vector_to_str(end_state_2_->d()) << "," << matrix_to_str(end_state_2_->r()) << "," << cutoff_ << "," << score_;
            out << std::endl;
        }
        
        if(record_state_) {
            int last = 0;
            for(int a = 2; a < sampler_.mst()->size(); a++) {
                
                if(sampler_.mst()->get_node(a)->data()->cur_state->size() > 4) {
                    last = a;
                    break;
                }
                out_state << sampler_.mst()->get_node(a)->data()->cur_state->name() << "|";
                out_state << sampler_.mst()->get_node(a)->data()->cur_state->end_ids()[0] << ",";
            }
            
            for(int a = last+1; a < sampler_.mst()->size(); a++) {
                if(sampler_.mst()->get_node(a)->data()->cur_state->size() > 4) {
                    last = a;
                    break;
                }
                out_state << sampler_.mst()->get_node(a)->data()->cur_state->name() <<  "|";
                out_state << sampler_.mst()->get_node(a)->data()->cur_state->end_ids()[0] << ",";
            }

            out_state << score_ << std::endl;
        }
        
        if(record_all_) {
            int last = 0;
            int first = 1;

            for(int a = 2; a < sampler_.mst()->size(); a++) {
                if(sampler_.mst()->get_node(a)->data()->cur_state->size() > 4) {
                    last = a;
                    break;
                }
                
                if(first) {
                    out_all << vector_to_str(sampler_.mst()->get_node(a)->data()->cur_state->end_states()[0]->d()) <<  ",";
                    out_all << matrix_to_str(sampler_.mst()->get_node(a)->data()->cur_state->end_states()[0]->r()) <<  ",";
                }
                
                out_all << vector_to_str(sampler_.mst()->get_node(a)->data()->cur_state->end_states()[1]->d()) <<  ",";
                out_all << matrix_to_str(sampler_.mst()->get_node(a)->data()->cur_state->end_states()[1]->r()) <<  ",";
                first = 0;
            }
            first = 1;
            
            for(int a = last+1; a < sampler_.mst()->size(); a++) {
                if(sampler_.mst()->get_node(a)->data()->cur_state->size() > 4) {
                    last = a;
                    break;
                }
                
                if(first) {
                    out_all << vector_to_str(sampler_.mst()->get_node(a)->data()->cur_state->end_states()[0]->d()) <<  ",";
                    out_all << matrix_to_str(sampler_.mst()->get_node(a)->data()->cur_state->end_states()[0]->r()) <<  ",";
                }
                
                out_all << vector_to_str(sampler_.mst()->get_node(a)->data()->cur_state->end_states()[1]->d()) <<  ",";
                out_all << matrix_to_str(sampler_.mst()->get_node(a)->data()->cur_state->end_states()[1]->r()) <<  ",";
                first = 0;
            }
            
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

    if(bound_pdbs_) {
        out_motifs.close();
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
    options_.add_option("steric_radius", 2.2f, OptionType::FLOAT);
    options_.add_option("all_pdbs", false, OptionType::BOOL);
    options_.add_option("bound_pdbs", false, OptionType::BOOL);
    options_.add_option("unbound_pdbs", false, OptionType::BOOL);
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
    temperature_    = options_.get_float("temperature");
    steps_          = options_.get_int("steps");
    cutoff_         = options_.get_float("cutoff");
    record_         = options_.get_bool("record");
    record_all_     = options_.get_bool("record_all");
    record_all_file_= options_.get_string("record_all_file");
    record_file_    = options_.get_string("record_file");
    record_state_   = options_.get_bool("record_state");
    steric_radius_  = options_.get_float("steric_radius");
    all_pdbs_       = options_.get_bool("all_pdbs");
    bound_pdbs_     = options_.get_bool("bound_pdbs");
    unbound_pdbs_   = options_.get_bool("unbound_pdbs");

}
