//
//  sequence_optimizer_3d.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/30/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef sequence_optimizer_3d_hpp
#define sequence_optimizer_3d_hpp

#include <stdio.h>

#include "base/types.h"
#include "base/option.h"
#include "util/random_number_generator.h"
#include "eternabot/scorer.h"
#include "motif_data_structure/motif_tree.h"
#include "motif_data_structure/motif_graph.h"
#include "motif_data_structure/motif_state_tree.h"
#include "motif_data_structure/motif_state_graph.hpp"


class SequenceOptimizerScorer {
public:
    SequenceOptimizerScorer(
            bool target_an_aligned_end) :
            target_an_aligned_end_(target_an_aligned_end) {}

public:
    virtual
    float
    score(
            motif_data_structure::MotifStateGraphOP const &) = 0;

public:
    float
    motif_state_diff(
            structure::BasepairStateOP const & end1,
            structure::BasepairStateOP const & end2) {
        auto diff = end1->d().distance(end2->d());
        if (target_an_aligned_end_) {
            diff += end1->r().difference(end2->r()) * 2;
        } else {
            end2->flip();
            diff += end1->r().difference(end2->r()) * 2;
            end2->flip();
        }
        return diff;
    }

protected:
    bool target_an_aligned_end_;

};

typedef std::shared_ptr<SequenceOptimizerScorer> SequenceOptimizerScorerOP;


class ExternalTargetScorer : public SequenceOptimizerScorer {
public:
    ExternalTargetScorer(
            structure::BasepairStateOP const & target,
            int ni,
            int ei,
            bool target_an_aligned_end) : SequenceOptimizerScorer(target_an_aligned_end),
            target_(target),
            ni_(ni),
            ei_(ei) {}

public:
    float
    score(
            motif_data_structure::MotifStateGraphOP const & msg) {
        state_ = msg->get_node(ni_)->data()->get_end_state(ei_);
        return motif_state_diff(state_, target_);
    }

private:
    structure::BasepairStateOP target_, state_;
    int ni_, ei_;

};


class InternalTargetScorer : public SequenceOptimizerScorer {
public:
    InternalTargetScorer(
            int ni1,
            int ei1,
            int ni2,
            int ei2,
            bool target_an_aligned_end) : SequenceOptimizerScorer(target_an_aligned_end),
            ni1_(ni1),
            ei1_(ei1),
            ni2_(ni2),
            ei2_(ei2) {}

public:
    float
    score(
            motif_data_structure::MotifStateGraphOP const & msg) {
        state1_ = msg->get_node(ni1_)->data()->get_end_state(ei1_);
        state2_ = msg->get_node(ni2_)->data()->get_end_state(ei2_);
        //return state1_->diff(state2_);
        return motif_state_diff(state1_, state2_);

    }

private:
    int ni1_, ni2_;
    int ei1_, ei2_;
    structure::BasepairStateOP state1_, state2_;
};

class MultiTargetScorer : public SequenceOptimizerScorer {
public:
    MultiTargetScorer(
            std::vector<SequenceOptimizerScorerOP> sub_scorers) :
            sub_scorers_(sub_scorers),
            SequenceOptimizerScorer(false) {}

public:
    float
    score(
            motif_data_structure::MotifStateGraphOP const & msg) {
        score_ = 0;
        for (auto const & sub_scorer : sub_scorers_) {
            score_ += sub_scorer->score(msg);
        }
        return score_;
    }

private:
    std::vector<SequenceOptimizerScorerOP> sub_scorers_;
    float score_;
};


class SequenceOptimizer3D {
public:

    SequenceOptimizer3D();

    ~SequenceOptimizer3D() {}

public: //setup

    void
    set_scorer(SequenceOptimizerScorerOP const & scorer) {
        scorer_ = scorer;
    }

private:

    struct OptimizedSequence {
        String sequence;
        float dist_score, eterna_score;
    };

    typedef std::shared_ptr<OptimizedSequence> OptimizedSequenceOP;
    typedef std::vector<OptimizedSequenceOP> OptimizedSequenceOPs;

    struct DesignableBP {
        inline
        DesignableBP(
                secondary_structure::BasepairOP const & nbp) :
                bp(nbp),
                last_state(Strings{"", ""}),
                m_id_bot(nullptr),
                m_id_top(nullptr) {}

        void
        update_state(
                Strings const & bp_name) {
            last_state[0] = bp->res1()->name();
            last_state[1] = bp->res2()->name();
            bp->res1()->name(bp_name[0]);
            bp->res2()->name(bp_name[1]);
        }

        void
        revert_state() {
            bp->res1()->name(last_state[0]);
            bp->res2()->name(last_state[1]);
        }


        secondary_structure::BasepairOP bp;
        Strings last_state;
        std::shared_ptr<util::Uuid> m_id_bot, m_id_top;

    };

    typedef std::shared_ptr<DesignableBP> DesignableBPOP;
    typedef std::vector<DesignableBPOP> DesignableBPOPs;

public:

    inline
    OptimizedSequenceOPs
    get_optimized_sequences(
            motif_data_structure::MotifGraphOP const & mg,
            SequenceOptimizerScorerOP const & scorer) {
        set_scorer(scorer);
        return get_optimized_sequences(mg);
    }

    OptimizedSequenceOPs
    get_optimized_sequences(
            motif_data_structure::MotifGraphOP const &);

    inline
    motif_data_structure::MotifGraphOP
    get_optimized_mg(
            motif_data_structure::MotifGraphOP const & mg,
            SequenceOptimizerScorerOP const & scorer) {
        set_scorer(scorer);
        return get_optimized_mg(mg);
    }

    motif_data_structure::MotifGraphOP
    get_optimized_mg(
            motif_data_structure::MotifGraphOP const &);

private:
    void
    _update_designable_bp(
            DesignableBPOP const &,
            motif_data_structure::MotifStateGraphOP &,
            secondary_structure::PoseOP &);

    String
    _validate_sequence(
            motif_data_structure::MotifStateGraphOP const &,
            secondary_structure::PoseOP const &);

    DesignableBPOPs
    _get_designable_bps(
            secondary_structure::PoseOP &);

    void
    _initiate_sequence_in_msg(
            motif_data_structure::MotifStateGraphOP &,
            secondary_structure::PoseOP const &);

    int
    convert_char_to_res_code(
            char c) {
        if (c == 'A') { return 0; }
        else if (c == 'C') { return 1; }
        else if (c == 'G') { return 2; }
        else if (c == 'U') { return 3; }
        else if (c == 'T') { return 3; }
        else if (c == 'N') { return -1; }
        else {
            throw secondary_structure::Exception("incorrect character for secondary string");
        }
    }

    void
    find_seq_violations(
            secondary_structure::PoseOP,
            Ints &);

    int
    find_gc_helix_stretches(
            secondary_structure::PoseOP);

    bool
    new_seq_violations() {
        for (int i = 0; i < current_violations_.size(); i++) {
            if (current_violations_[i] != next_violations_[i]) { return true; }
        }

        if (current_gc_stretches_ < next_gc_stretches_) { return true; }

        return false;
    }

public: //option wrappers

    inline
    base::Options &
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

    template<typename T>
    void
    set_option_value(
            String const & name,
            T const & val) {
        options_.set_value(name, val);
        update_var_options();
    }

    void
    update_var_options();

protected:
    void
    setup_options();

private:
    base::Options options_;
    eternabot::Scorer eterna_scorer_;
    util::RandomNumberGenerator rng_;
    SequenceOptimizerScorerOP scorer_;
    std::vector<Strings> possible_bps_;
    // option vars
    int solutions_;
    int steps_;
    float cutoff_, eterna_cutoff_;
    bool verbose_, return_lowest_;
    Strings disallowed_sequences_;
    std::vector<Ints> disallowed_res_types_sequences_;
    Ints current_violations_;
    Ints next_violations_;
    int current_gc_stretches_, next_gc_stretches_;


};

typedef std::shared_ptr<SequenceOptimizer3D> SequenceOptimizer3DOP;

#endif /* sequence_optimizer_3d_hpp */
