//
//  sequence_designer.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include <base/log.hpp>
#include <secondary_structure/sequence_tools.h>
#include <eternabot/sequence_designer.h>

namespace eternabot {


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// setup functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SequenceDesigner::SequenceDesigner():
        scorer_(Scorer()),
        results_(SequenceDesignerResultOPs()),
        rng_(util::RandomNumberGenerator()) {

    possible_bps_ = std::vector<Strings>({{"A", "U"}, {"U", "A"}, {"G", "C"}, {"C", "G"}});


    // generate constraints
    auto disallowed_sequences = Strings{"AAAA", "CCCC", "GGGG", "UUUU"};
    for(auto const & seq : disallowed_sequences) { seq_constraints_.add_disallowed_sequence(seq); }
    seq_constraints_.add_gc_helix_stretch_limit(3);
    current_violations_ = Ints(seq_constraints_.num_constraints());
    next_violations_    = Ints(seq_constraints_.num_constraints());

    setup_options();
    temperature_ = 4.0f;
    mc_ = util::MonteCarlo(temperature_);
}


void
SequenceDesigner::setup() {
    /*if((int)results_.size() != steps_) { results_.resize(steps_); }
    for(auto & r : results_) {
        r = std::make_shared<SequenceDesignerResult>();
    }*/



}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// option functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
SequenceDesigner::setup_options() {
    options_.add_option("designs", 1, base::OptionType::INT);
    options_.add_option("steps", 1000, base::OptionType::INT);
    options_.lock_option_adding();
    update_var_options();
    
}
    
void
SequenceDesigner::update_var_options() {
    designs_ = options_.get_int("designs");
    steps_   = options_.get_int("steps");
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
SequenceDesignerResultOPs const &
SequenceDesigner::design(
        secondary_structure::PoseOP const & p) {
    results_ = SequenceDesignerResultOPs();
    _find_designable_bps(p);       // find all basepairs with N-N residues
    _generate_inital_sequence(p);  // generate initial sequence, get rid as many violations as possible
    scorer_.setup(p);

    auto current_score = scorer_.score_secondary_structure(p);
    auto current_sequence = p->sequence();
    auto best_score = current_score;
    auto best_sequence = p->sequence();
    auto next_score = 0.0f;
    auto last_pair = Strings{"", ""};

    for(int i = 0; i < steps_; i++) {
        auto pos = rng_.randrange(designable_bps_.size());
        auto & pair = _get_random_pair();
        last_pair[0] = designable_bps_[pos]->res1()->name();
        last_pair[1] = designable_bps_[pos]->res2()->name();
        _set_bp_sequence(pair, designable_bps_[pos]);
        next_violations_ = seq_constraints_.violations(p);

        if(_new_sequence_violations()) {
            _set_bp_sequence(last_pair, designable_bps_[pos]);
            continue;
        }

        next_score = scorer_.score_secondary_structure(p);
        if(mc_.accept(current_score, next_score)) {
            current_score = next_score;
        }
        else {
            _set_bp_sequence(last_pair, designable_bps_[pos]);
            continue;
        }

        if(current_score > best_score) {
            best_score = current_score;
            best_sequence = p->sequence();
        }
    }

    results_.push_back(std::make_shared<SequenceDesignerResult>(best_sequence, best_score));
    return results_;
    
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Strings const &
SequenceDesigner::_get_random_pair() {
    return possible_bps_[rng_.randrange(possible_bps_.size())];
}

void
SequenceDesigner::_set_bp_sequence(
        Strings const & pair,
        secondary_structure::BasepairOP bp) {
    bp->res1()->name(pair[0]);
    bp->res2()->name(pair[1]);
}


void
SequenceDesigner::_find_designable_bps(
        secondary_structure::PoseOP p) {
    for(auto const & bp : p->basepairs()) {
        if(bp->res1()->res_type() == secondary_structure::ResType::NONE &&
           bp->res2()->res_type() == secondary_structure::ResType::NONE ) {
            designable_bps_.push_back(bp);
        }
    }

    if(designable_bps_.size() == 0) {
        LOG_ERROR << "no basepairs are designable! "; exit(0);
    }

}


bool
SequenceDesigner::_new_sequence_violations() {
    for(int i = 0; i < current_violations_.size(); i++) {
        if(current_violations_[i] != next_violations_[i]) { return true; }
    }
    return false;
}

void
SequenceDesigner::_generate_inital_sequence(
        secondary_structure::PoseOP p) {

    auto count = 0;
    current_violations_ = seq_constraints_.violations(p);

    // inital sequence filling
    for (auto & bp : designable_bps_) {
        auto & p = _get_random_pair();
        _set_bp_sequence(p, bp);
    }

    next_violations_ = seq_constraints_.violations(p);
    while(_new_sequence_violations()) {
        count += 1;
        auto pos = rng_.randrange(designable_bps_.size());
        auto & pair = _get_random_pair();
        _set_bp_sequence(pair, designable_bps_[pos]);
        if(count > 1000000) {
            LOG_WARNING << "cannot find initial sequence that does not have sequence violations! ";
            current_violations_ = next_violations_;
            break;
        }
        next_violations_ = seq_constraints_.violations(p);
    }

}

}























