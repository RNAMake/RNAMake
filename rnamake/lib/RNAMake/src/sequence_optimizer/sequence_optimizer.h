//
//  sequence_optimizer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/11/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sequence_optimizer__
#define __RNAMake__sequence_optimizer__

#include <stdio.h>
#include "base/option.h"
#include "eternabot/sequence_designer.h"
#include "motif_data_structures/motif_graph.h"
#include "motif_data_structures/motif_tree.h"

struct SequenceOptimizerResult {
    inline
    SequenceOptimizerResult(
        MotifTreeOP const & nmt,
        float nscore):
    motif_tree(nmt),
    score(nscore)
    {}
    
    MotifTreeOP motif_tree;
    float score;
};

struct OptimizedSequence {
    inline
    OptimizedSequence(
        String const & nsequence,
        float const & nclose_distance,
        float const & neternabot_score):
    sequence(nsequence),
    close_distance(nclose_distance),
    eternabot_score(neternabot_score)
    {}
    
    String sequence;
    float close_distance, eternabot_score;
};

typedef std::shared_ptr<SequenceOptimizerResult> SequenceOptimizerResultOP;
typedef std::shared_ptr<OptimizedSequence>       OptimizedSequenceOP;
typedef std::vector<OptimizedSequenceOP>         OptimizedSequenceOPs;


class SequenceOptimizer {
public:
    SequenceOptimizer():
    options_(Options()){ setup_options(); }
    
    ~SequenceOptimizer() {}
    
public:

    
    OptimizedSequenceOPs
    get_optimized_sequences(
        MotifGraphOP & mg,
        Uuid const & uuid_1,
        Uuid const & uuid_2,
        int end_i,
        int end_j);
    
    /*OptimizedSequenceOPs
    get_optimized_sequences_section(
        MotifGraphOP & mg,
        MotifGraphOP & section,
        Uuid const & uuid_2,
        int end_i,
        int end_j);*/
    
    
    String
    get_final_sequence(
        String const &,
        String const &);
    
public: //option wrappers
    
    inline
    Options &
    options() { return options_; }
    
    inline
    float
    get_int_option(String const & name) { return options_.get_int(name); }
    
    inline
    float
    get_float_option(String const & name) { return options_.get_float(name); }
    
    inline
    String
    get_string_option(String const & name) { return options_.get_string(name); }
    
    inline
    bool
    get_bool_option(String const & name) { return options_.get_bool(name); }
    
    inline
    bool
    has_option(String const & name) { return options_.has_option(name); }
    
    template<typename T>
    void
    set_option_value(
        String const & name,
        T const & val) {
        options_.set_value(name, val);
        update_var_options();
    }
    
protected:
    
    void
    setup_options();
    
    void
    update_var_options();
    
    
private:
    Options options_;
    MotifTreeOP mt_;
    eternabot::SequenceDesigner designer_;
    eternabot::SequenceDesignerResultOPs designer_results_;
    //options
    bool sub_sequence_;
    int start_;
    int end_;
    int solutions_;
};

#endif /* defined(__RNAMake__sequence_optimizer__) */










