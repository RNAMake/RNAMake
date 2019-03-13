//
//  mini_ttr.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/16/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__mini_ttr__
#define __RNAMake__mini_ttr__

#include <stdio.h>

//RNAMake Headers
#include "base/option.h"
#include "base/cl_option.h"
#include "motif_data_structure/motif_graph.h"
#include "motif_search/motif_state_search.h"
#include "resources/resource_manager.h"
#include "sequence_optimization/sequence_optimizer.h"

base::CommandLineOptions
parse_command_line(
    int argc,
    const char ** argv);

class MiniTTR {
public:
    MiniTTR():
    search_(motif_search::MotifStateSearch()),
    optimizer_(SequenceOptimizer()),
    options_(Options("MiniTTROptions"))
    {}
    
    ~MiniTTR() {}

public:
    
    virtual
    void
    setup(base::CommandLineOptions const & opts) {
        setup_options();
        search_.set_option_value("max_node_level", 10);
        search_.set_option_value("min_node_level", 0);
        search_.set_option_value("max_solutions", 1000);
        search_.set_option_value("accept_score", 10);
        search_.set_option_value("min_ss_score", 0);
        search_.set_option_value("max_size", 100);
        
        for(auto const & opt : opts) {
            if(! search_.has_option(opt->name())) { continue; }
            if(! opts.is_filled(opt->name())) { continue; }
            if     (opt->type() == base::OptionType::INT) {
                search_.set_option_value(opt->name(), opt->get_int());
            }
            else if(opt->type() == base::OptionType::FLOAT) {
                search_.set_option_value(opt->name(), opt->get_float());
            }
            else if(opt->type() == base::OptionType::BOOL) {
                search_.set_option_value(opt->name(), opt->get_bool());
            }
            else if(opt->type() == base::OptionType::STRING) {
                search_.set_option_value(opt->name(), opt->get_string());
            }
        }
        
        options_.set_value("test_run", opts.get_bool("test_run"));
        options_.set_value("out", opts.get_string("out"));
        update_var_options();
        
        /*for(auto const & opt: opts) {
         if(! options_.has_option(opt->name())) { continue; }
         if(! opts.is_filled(opt->name())) { continue; }
         }*/
        
        mg_ = motif_data_structure::MotifGraph();
        /*mg_.add_motif("HELIX.IDEAL.2");
        mg_.add_motif("GAAA_tetraloop", "A229-A245");
        mg_.add_motif("HELIX.IDEAL.3", -1, "A149-A154");
        mg_.add_motif("HELIX.IDEAL.3", 1, "A222-A251");*/
        mg_.write_pdbs();
        
        mg_.set_option_value("sterics", false);
        
    }
    
    virtual
    void
    run();
    
protected:
    Options options_;
    motif_search::MotifStateSearch search_;
    SequenceOptimizer optimizer_;
    motif_data_structure::MotifGraph mg_;
    bool test_run_, opt_seq_;

private:
    void
    optimize_sequence(motif_data_structure::MotifGraph &);
    
    virtual
    void
    setup_options() {
        options_.add_option("test_run", false, base::OptionType::BOOL);
        options_.add_option("out", String("solutions.top"), base::OptionType::STRING);
        options_.add_option("opt_seq", true, base::OptionType::BOOL);
        options_.lock_option_adding();
        update_var_options();
    }
    
    void
    update_var_options() {
        test_run_ = options_.get_bool("test_run");
        opt_seq_  = options_.get_bool("opt_seq");
    }

    
};


class MiniTTRPathFollow : public MiniTTR {
public:
    
    MiniTTRPathFollow() : MiniTTR() {
    }

    ~MiniTTRPathFollow() {}
    
public:
    void
    setup(base::CommandLineOptions const & opts) {
        setup_options();
        search_.set_option_value("max_node_level", 400);
        search_.set_option_value("min_node_level", 0);
        search_.set_option_value("max_solutions", 100000000);
        search_.set_option_value("accept_score", 15);
        search_.set_option_value("max_size", 10000);
        
        for(auto const & opt : opts) {
            if(! search_.has_option(opt->name())) { continue; }
            if(! opts.is_filled(opt->name())) { continue; }
            if     (opt->type() == base::OptionType::INT) {
                search_.set_option_value(opt->name(), opt->get_int());
            }
            else if(opt->type() == base::OptionType::FLOAT) {
                search_.set_option_value(opt->name(), opt->get_float());
            }
            else if(opt->type() == base::OptionType::BOOL) {
                search_.set_option_value(opt->name(), opt->get_bool());
            }
            else if(opt->type() == base::OptionType::STRING) {
                search_.set_option_value(opt->name(), opt->get_string());
            }
        }
        
        options_.set_value("path", opts.get_string("path"));
        options_.set_value("test_run", opts.get_bool("test_run"));
        std::cout << "made it" << std::endl;
        update_var_options();

        /*for(auto const & opt: opts) {
            if(! options_.has_option(opt->name())) { continue; }
            if(! opts.is_filled(opt->name())) { continue; }
        }*/
        
        mg_ = motif_data_structure::MotifGraph();
        mg_.set_option_value("sterics", false);
        /*mg_.add_motif("GAAA_tetraloop", "A229-A245");
        mg_.add_motif("HELIX.IDEAL.6", -1, "A149-A154");*/
    }
    
    void
    run();
    
private:
    
    void
    setup_options() {
        options_.add_option("path", String(""), base::OptionType::STRING);
        options_.add_option("test_run", false, base::OptionType::BOOL);
        options_.lock_option_adding();
        update_var_options();
    }
    
    void
    update_var_options() {
        test_run_ = options_.get_bool("test_run");
    }
    
};



#endif /* defined(__RNAMake__mini_ttr__) */
