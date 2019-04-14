//
//  sequence_designer.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include <base/log.h>
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

    // generate disallowed sequences
    auto disallowed_sequences = Strings{"AAAA", "CCCC", "GGGG", "UUUU"};
    for(auto const & seq : disallowed_sequences) {
        auto disallowed_types = secondary_structure::ResTypes();
        secondary_structure::get_res_types_from_sequence(seq, disallowed_types);
        disallowed_res_type_arrays_.push_back(disallowed_types);
    }

    current_violations_ = Ints(disallowed_sequences.size());
    next_violations_ = Ints(disallowed_sequences.size());

    setup_options();
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
        secondary_structure::PoseOP const & p_org) {

    auto p = std::make_shared<secondary_structure::Pose>(*p_org);
    // find all basepairs with N-N types
    _find_designable_bps(p);
    // generate initial sequence, get rid as many violations as possible
    _generate_inital_sequence(p);

    /* scorer_.setup(p);

     designable_bps_ = secondary_structure::BasepairOPs();
     temperature_ = 4.0f;

     exit(0);


     for(auto const & bp : p->basepairs()) {
         if     (bp->res1()->res_type() == -1 && bp->res2()->res_type() == -1) {
             designable_bps_.push_back(bp); }
         else if(bp->res1()->res_type() == -1 || bp->res2()->res_type() == -1) {
             throw std::runtime_error("only one residue in a basepair is marked as N, this is not allowed");
         }
     }

     auto pair = Strings();
     for(auto & bp : designable_bps_) {
         pair = pairs[rng_.randrange((int)pairs.size())];
         bp->res1()->name(pair[0]);
         bp->res2()->name(pair[1]);
     }

     auto res = p->residues();
     for(auto & bp : designable_bps_) {
         if(_row_of_gc_bps(bp->res1(), res) == 1 || _row_of_gc_bps(bp->res2(), res) == 1) {
             if (bp->res1()->name() == "G" && bp->res2()->name() == "C") {
                 bp->res1()->name("A"); bp->res2()->name("U");
             }
             else if(bp->res1()->name() == "C" && bp->res2()->name() == "G") {
                 bp->res1()->name("U"); bp->res2()->name("A");

             }
         }

     }

     std::cout << p->sequence() << std::endl;

     auto last_score = scorer_.score_secondary_structure(p);
     auto current_score = 0.0f;
     String last_sequence = p->sequence();
     int pos = 0;
     int pair_pos = 0;
     int j = 0;
     int interval =  steps_ / 10 ;
     float prob, diceroll;
     int fail = 0;

     String best_sequence = p->sequence();
     float best_score = -100000;

     Strings last_pair{"", ""};
     for(int i = 0; i < steps_; i++) {
         if(i > 0 && i % interval == 0) {
             temperature_ = temperature_ * 0.8;
             p->replace_sequence(best_sequence);
             last_sequence = best_sequence;
             last_score = best_score;
         }

         pos = rng_.randrange((int)designable_bps_.size());

         fail = 1;
         if((designable_bps_[pos]->res1()->name() == "G" && designable_bps_[pos]->res2()->name() == "C") ||
            (designable_bps_[pos]->res1()->name() == "C" && designable_bps_[pos]->res2()->name() == "G")) {

             if(_row_of_gc_bps(designable_bps_[pos]->res1(), res) == 1 ||
                _row_of_gc_bps(designable_bps_[pos]->res2(), res) == 1) {
                 diceroll = rng_.rand();
                 if(diceroll < 0.5) { pair = pairs[0]; }
                 else               { pair = pairs[1]; }
                 fail = 0;
             }

         }

         while(fail) {
             pair_pos = rng_.randrange((int)pairs.size());
             if(pair_pos == 4) { pair_pos = 0; }
             pair = pairs[pair_pos];
             if(pair[0] == designable_bps_[pos]->res1()->name() &&
                pair[1] == designable_bps_[pos]->res2()->name()) {
                 continue;
             }
             break;
         }
         last_pair[0] = designable_bps_[pos]->res1()->name();
         last_pair[1] = designable_bps_[pos]->res2()->name();
         designable_bps_[pos]->res1()->name(pair[0]);
         designable_bps_[pos]->res2()->name(pair[1]);
         current_score = scorer_.score_secondary_structure(p);

         results_[j]->sequence = p->sequence();
         results_[j]->score = current_score;
         j++;

         if(current_score > last_score) {
             last_score = current_score;
             last_sequence = results_[j-1]->sequence;

             if(current_score > best_score) {
                 best_score = current_score;
                 best_sequence = results_[j-1]->sequence;
             }

             continue;
         }

         prob = expf((current_score - last_score) / temperature_);
         diceroll = rng_.rand();
         if(diceroll < prob) {
             last_score = current_score;
             last_sequence = results_[j-1]->sequence;
             continue;
         }

         designable_bps_[pos]->res1()->name(last_pair[0]);
         designable_bps_[pos]->res2()->name(last_pair[1]);

     }

     std::sort(results_.begin(), results_.end(),
               sequence_designer_result_less_than_key());
     */
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
    bp->res1()->name(pair[1]);
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

void
SequenceDesigner::_generate_inital_sequence(
        secondary_structure::PoseOP p) {

    

    for(auto & bp : designable_bps_) {

    }

}

}























