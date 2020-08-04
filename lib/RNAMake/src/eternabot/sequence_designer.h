//
//  sequence_designer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sequence_designer__
#define __RNAMake__sequence_designer__

#include <stdio.h>

#include <base/option.h>
#include <util/random_number_generator.h>
#include <util/monte_carlo.h>
#include <secondary_structure/sequence_constraint.h>
#include <eternabot/scorer.h>

namespace eternabot {

struct SequenceDesignerResult {
    inline
    SequenceDesignerResult(
            String const & n_sequence,
            float n_score):
            sequence(n_sequence),
            score(n_score) {}

    String sequence;
    float score;
};
    
typedef std::shared_ptr<SequenceDesignerResult> SequenceDesignerResultOP;
typedef std::vector<SequenceDesignerResultOP> SequenceDesignerResultOPs;

struct sequence_designer_result_less_than_key {
    inline bool operator() (
        SequenceDesignerResultOP const & r1,
        SequenceDesignerResultOP const & r2) {
    
        return (r1->score > r2->score);
    }
};

    
class SequenceDesigner {
public:
    SequenceDesigner();
    
    ~SequenceDesigner() {}

public:
    
    void
    setup();
    
    SequenceDesignerResultOPs const &
    design(secondary_structure::PoseOP const &);
    
    
public: //option wrappers

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
    Strings const &
    _get_random_pair();

    void
    _set_bp_sequence(
            Strings const &,
            secondary_structure::BasepairOP);

    void
    _find_designable_bps(
            secondary_structure::PoseOP);

    void
    _generate_inital_sequence(
            secondary_structure::PoseOP);

    bool
    _new_sequence_violations();

private:
    struct Parameters {

    };


private:
    std::vector<Strings> possible_bps_;
    base::Options options_;
    util::RandomNumberGenerator rng_;
    util::MonteCarlo mc_;
    secondary_structure::BasepairOPs designable_bps_;
    Scorer scorer_;
    SequenceDesignerResultOPs results_;

    // tracking sequence constraints
    Ints current_violations_, next_violations_;
    secondary_structure::SequenceConstraints seq_constraints_;

    int designs_, steps_;
    float temperature_;
    
    
};
    
}

#endif /* defined(__RNAMake__sequence_designer__) */
