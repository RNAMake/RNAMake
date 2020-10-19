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
            float n_score,
            float n_bp_diff_score):
            sequence(n_sequence),
            score(n_score),
            bp_diff_score(n_bp_diff_score){}

    String sequence;
    float score;
    float bp_diff_score;
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

class MonteCarloMove {
public:
    MonteCarloMove() = default;

    virtual
    ~MonteCarloMove() = default;

    virtual
    MonteCarloMove *
    clone() const = 0;

public:
    virtual
    int
    move(
            secondary_structure::PoseOP) = 0;

    virtual
    void
    undo(
            secondary_structure::PoseOP) = 0;
};

class MutateBPMove : public MonteCarloMove {
public:
    MutateBPMove(
            secondary_structure::BasepairOPs const & designable_bps,
            std::vector<secondary_structure::ResTypes> const & possible_rt_types,
            std::map<int, secondary_structure::ResType> const & res_type_constraints,
            std::map<secondary_structure::BasepairOP, int> const & closing_bps):
            designable_bps_(designable_bps),
            possible_rt_types_(possible_rt_types),
            res_type_constraints_(res_type_constraints),
            closing_bps_(closing_bps),
            rng_(util::RandomNumberGenerator()) {
        last_res_types_ = secondary_structure::ResTypes(2);
    }

    ~MutateBPMove() = default;

    MonteCarloMove *
    clone() const { return new MutateBPMove(*this); }

public:
    int
    move(
            secondary_structure::PoseOP p) {

        if(designable_bps_.size() == 0) { return 0; }

        current_ = designable_bps_[rng_.randrange((int)designable_bps_.size())];
        last_res_types_[0] = current_->res1()->res_type();
        last_res_types_[1] = current_->res2()->res_type();
        int count = 0;

        auto org_res_type_1 = res_type_constraints_[current_->res1()->num()];
        auto org_res_type_2 = res_type_constraints_[current_->res2()->num()];

        while(true) {
            count += 1;
            if(count > 1000) { break; }
            if(rng_.randrange(1000) > 200) {
                current_res_types_ = possible_rt_types_[rng_.randrange(4)];
            }
            else {
                current_res_types_ = possible_rt_types_[rng_.randrange(6)];

            }

            if(closing_bps_.find(current_) != closing_bps_.end()) {
                _get_random_res_type_pair_gc_cap(current_res_types_);
            }

            //std::cout << current_res_types_[0] << " " << current_res_types[1];
            if(!secondary_structure::does_restype_satisfy_constraint(current_res_types_[0], org_res_type_1) ||
               !secondary_structure::does_restype_satisfy_constraint(current_res_types_[1], org_res_type_2)) {
                continue;
            }

            // got the same pair again, should not count that as a successful move - JDY
            if(current_res_types_[0] == current_->res1()->res_type()) { return 0;}

            current_->res1()->res_type(current_res_types_[0]);
            current_->res2()->res_type(current_res_types_[1]);

            return 1;

        }

        return 0;
    }

    void
    undo(
            secondary_structure::PoseOP p) {
        current_->res1()->res_type(last_res_types_[0]);
        current_->res2()->res_type(last_res_types_[1]);
    }

private:
    void
    _get_random_res_type_pair_gc_cap(
            secondary_structure::ResTypes & pair) {
        // biased base pair selection for caped ends, 80% chance to be GC/CG over AU/UA
        // will do GCs
        auto rand = rng_.randrange(1000);
        if(rand > 200) {
            // selecting GC
            if(rng_.randrange(1000) > 500) {
                pair[0] = secondary_structure::ResType::G;
                pair[1] = secondary_structure::ResType::C;
            }
                // selecting CG
            else {
                pair[0] = secondary_structure::ResType::C;
                pair[1] = secondary_structure::ResType::G;
            }
        }

        else if(rand > 50){
            // selecting AU
            if(rng_.randrange(1000) > 500) {
                pair[0] = secondary_structure::ResType::A;
                pair[1] = secondary_structure::ResType::U;
            }
                // selecting UA
            else {
                pair[0] = secondary_structure::ResType::U;
                pair[1] = secondary_structure::ResType::A;
            }
        }
        else {
            // selecting AU
            if(rng_.randrange(1000) > 500) {
                pair[0] = secondary_structure::ResType::G;
                pair[1] = secondary_structure::ResType::U;
            }
                // selecting UA
            else {
                pair[0] = secondary_structure::ResType::U;
                pair[1] = secondary_structure::ResType::G;
            }

        }
    }


private:
    secondary_structure::BasepairOPs designable_bps_;
    secondary_structure::BasepairOP current_;
    std::vector<secondary_structure::ResTypes> possible_rt_types_;
    std::map<int, secondary_structure::ResType> res_type_constraints_;
    std::map<secondary_structure::BasepairOP, int> closing_bps_;
    util::RandomNumberGenerator rng_;
    secondary_structure::ResTypes last_res_types_, current_res_types_;

};

class MutateUnpairedResMove : public MonteCarloMove {
public:
    MutateUnpairedResMove(
            secondary_structure::ResidueOPs designable_unpaired_res,
            std::map<int, secondary_structure::ResType> const & res_type_constraints):
            designable_unpaired_res_(designable_unpaired_res),
            res_type_constraints_(res_type_constraints) {
        possible_res_types_ = secondary_structure::ResTypes {
            secondary_structure::ResType::A,
            secondary_structure::ResType::C,
            secondary_structure::ResType::G,
            secondary_structure::ResType::U
        };

    }

    ~MutateUnpairedResMove() override = default;

    MonteCarloMove *
    clone() const override { return new MutateUnpairedResMove(*this); }

public:
    int
    move(
            secondary_structure::PoseOP p) override {

        if(designable_unpaired_res_.empty()) { return 0; }
        current_ = designable_unpaired_res_[rng_.randrange((int)designable_unpaired_res_.size())];

        auto org_res_type = res_type_constraints_[current_->num()];

        last_res_type_ = current_->res_type();
        auto count = 0;
        while(true) {
            count += 1;
            if(count > 1000) { return 0; }
            current_res_type_ = possible_res_types_[rng_.randrange((int)possible_res_types_.size())];
            /*if(org_res_type == secondary_structure::ResType::N) {
                auto rand = rng_.randrange(1000);
                if(rand < 900) {
                    current_res_type_ = secondary_structure::ResType::A;
                }
            }*/

            if(!secondary_structure::does_restype_satisfy_constraint(current_res_type_, org_res_type)) {
                continue;
            }
            if(current_res_type_ == last_res_type_) { return 0; }
            current_->res_type(current_res_type_);
            break;
        }
        return 1;
    }

    void
    undo(
            secondary_structure::PoseOP p) override {
        current_->res_type(last_res_type_);

    }

private:
    secondary_structure::ResidueOP current_;
    secondary_structure::ResidueOPs designable_unpaired_res_;
    secondary_structure::ResTypes possible_res_types_;
    secondary_structure::ResType last_res_type_, current_res_type_;
    std::map<int, secondary_structure::ResType> res_type_constraints_;
    util::RandomNumberGenerator rng_;

};

/*class BoostLoopMove : public MonteCarloMove {


};*/


typedef std::shared_ptr<MonteCarloMove> MonteCarloMoveOP;
typedef std::vector<MonteCarloMoveOP> MonteCarloMoveOPs;

    
class SequenceDesigner {
public:
    SequenceDesigner();
    
    ~SequenceDesigner() = default;

public:
    
    void
    setup();
    
    SequenceDesignerResultOPs const &
    design(secondary_structure::PoseOP const &);

public:
    void
    set_previous_solutions(
            Strings const & solutions) {
        previous_solutions_ = solutions;
    }
    
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
    void
    _find_designable_bps(
            secondary_structure::PoseOP);

    void
    _generate_inital_sequence(
            secondary_structure::PoseOP);

    bool
    _new_sequence_violations(
            Ints const &,
            Ints const &);


private: // new and possibly badly made functions
    void
    _set_initial_helix_sequence(
            secondary_structure::MotifOP);

    void
    _get_random_res_type_pair_gc_cap(
            secondary_structure::ResTypes &);

    void
    _get_random_res_type_pair(
            secondary_structure::ResTypes &);

    float
    _optimize_substructure(
            secondary_structure::PoseOP,
            int);

    float
    _bp_list_diff(
            secondary_structure::PoseOP,
            std::vector<std::vector<int>> const &,
            size_t,
            FeaturesOP);

private:
    bool
    _assign_new_residue_restype(
            secondary_structure::ResidueOP);

private:
    struct Parameters {
        bool biased_gc_caps;
    };


private:
    secondary_structure::ResidueOPs designable_res_;
    secondary_structure::ResidueOPs designable_unpaired_res_;
    std::vector<Strings> possible_bps_;
    std::vector<secondary_structure::ResTypes> possible_rt_bps_;
    std::vector<secondary_structure::ResType> possible_res_types_;
    std::map<secondary_structure::BasepairOP, int> closing_bps_;
    base::Options options_;
    util::RandomNumberGenerator rng_;
    util::MonteCarlo mc_;
    std::map<int, secondary_structure::ResType> res_type_constraints_;
    secondary_structure::BasepairOPs designable_bps_;
    Scorer scorer_;
    SequenceDesignerResultOPs results_;

    // tracking sequence constraints
    Ints current_violations_, next_violations_;
    secondary_structure::SequenceConstraints seq_constraints_;

    // current solutions
    std::map<secondary_structure::MotifOP, secondary_structure::ResTypes> current_restypes_;

    // directly keep track of vienna pair map
    std::vector<std::vector<int>> pair_map_;
    size_t pair_map_entries_;
    vienna::Vienna v_;

    Strings previous_solutions_;


    int designs_, steps_;
    float temperature_;

    Parameters parameters_;
    
    
};
    
}

#endif /* defined(__RNAMake__sequence_designer__) */
