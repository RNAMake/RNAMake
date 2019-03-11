//
//  motif_state.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state__
#define __RNAMake__motif_state__

#include <stdio.h>

#include "base/types.h"
#include "math/xyz_vector.h"
#include "math/xyz_matrix.h"
#include "util/uuid.h"
#include "structure/basepair_state.fwd.h"
#include "structure/basepair_state.h"

/*
* Exception for MotifState
*/
class MotifStateException : public std::runtime_error {
public:
    /**
     * Standard constructor for MotifStateException
     * @param   message   Error message for motif state
     */
    MotifStateException(String const & message):
    std::runtime_error(message)
    {}
};


class MotifState {
public:
    inline
    MotifState(
        String const & name,
        Strings const & end_names,
        Strings const & end_ids,
        BasepairStateOPs const & end_states,
        math::Points const & beads,
        float const & score,
        int const & size,
        int const & block_end_add):
    name_(name),
    end_names_(end_names),
    end_ids_(end_ids),
    end_states_(end_states),
    beads_(beads),
    score_(score),
    size_(size),
    block_end_add_(block_end_add),
    uuid_(Uuid())
    {}
    
    inline
    MotifState(
        String const & name,
        Strings const & end_names,
        Strings const & end_ids,
        BasepairStateOPs const & end_states,
        math::Points const & beads,
        float const & score,
        int const & size,
        int const & block_end_add,
        Uuid const & uuid):
    name_(name),
    end_names_(end_names),
    end_ids_(end_ids),
    end_states_(end_states),
    beads_(beads),
    score_(score),
    size_(size),
    block_end_add_(block_end_add),
    uuid_(uuid)
    {}
    
    inline
    MotifState(
        MotifState const & ms):
    name_(ms.name_),
    end_names_(ms.end_names_),
    end_ids_(ms.end_ids_),
    end_states_(BasepairStateOPs(ms.end_states_.size())),
    beads_(ms.beads_),
    score_(ms.score_),
    size_(ms.size_),
    block_end_add_(ms.block_end_add_),
    uuid_(ms.uuid_) {
        int i = 0;
        for(auto const & bp_state : ms.end_states_) {
            end_states_[i] = std::make_shared<BasepairState>(*bp_state);
            i++;
        }
    }
    
    inline
    MotifState(
        String const & s) {
        
        auto spl = base::split_str_by_delimiter(s, "|");
        name_ = spl[0];
        score_ = std::stof(spl[1]);
        size_ = std::stoi(spl[2]);
        block_end_add_ = std::stoi(spl[3]);
        beads_ = math::vectors_from_str(spl[4]);
        end_names_ = base::split_str_by_delimiter(spl[5], ",");
        end_ids_ = base::split_str_by_delimiter(spl[6], ",");
        end_states_ = BasepairStateOPs();
        for(int i = 7; i < spl.size(); i++) {
            auto bp_state = std::make_shared<BasepairState>(spl[i]);
            end_states_.push_back(bp_state);
        }
        uuid_ = Uuid();

        
    }
    
    
    ~MotifState() {}

public: // non const methods
    inline
    void
    move(
            math::Point const & p) {
        for(auto & end_state : end_states_) {
            end_state->move(p);
        }

        for( auto & b : beads_) {
            b += p;
        }
    }

    inline
    void
    transform(
            math::Matrix const & r,
            math::Vector const & t,
            math::Point & dummy) {

        for(auto & end_state : end_states_) {
            end_state->transform(r, t, dummy);
        }

        for( auto & b : beads_) {
            math::dot_vector(r, b, dummy);
            dummy += t;
            b = dummy;
        }
    }

    inline
    void
    transform(
            math::Matrix const & r,
            math::Vector const & t) {
        auto dummy = math::Point();
        transform(r, t, dummy);
    }

public:
    
    
    String
    to_str();
    
    BasepairStateOP const &
    get_end_state(
        String const &);
    
    int
    get_end_index(
        BasepairStateOP const &);
    
    int
    get_end_index(
        String const &);
    
public: // getters
    
    inline
    String const &
    name() { return name_; }
    
    inline
    Strings const &
    end_names() { return end_names_; }
    
    inline
    Strings const &
    end_ids() { return end_ids_; }
    
    inline
    BasepairStateOPs  &
    end_states() { return end_states_; }
    
    inline
    math::Points const &
    beads() { return beads_; }
    
    inline
    float const &
    score() { return score_; }
    
    inline
    int const &
    size() { return size_; }
    
    inline
    int const &
    block_end_add() { return block_end_add_; }
    
    inline
    Uuid const &
    uuid() { return uuid_; }
    
public: //setters
  
    inline
    void
    uuid(Uuid const & uuid) { uuid_ = uuid; }
    
    
public:
    inline
    void
    update_end_state(
        int i,
        BasepairState const & new_state) {
        end_states_[i]->set(new_state);
    }
    
    inline
    void
    beads(
        math::Points const & beads) { beads_ = beads; }
    
    inline
    void
    new_uuids() {
        uuid_ = Uuid();
    }
    
    
    
private:
    String name_;
    Strings end_names_;
    Strings end_ids_;
    BasepairStateOPs end_states_;
    math::Points beads_;
    float score_;
    int size_;
    int block_end_add_;
    Uuid uuid_;
    
};

typedef std::shared_ptr<MotifState> MotifStateOP;
typedef std::vector<MotifStateOP>   MotifStateOPs;


inline
void
align_motif_state(
    BasepairStateOP const & ref_bp_state,
    MotifStateOP & org_state) {
    
    BasepairState bp_state, bp_state_final;
    
    ref_bp_state->get_transforming_r_and_t(*org_state->end_states()[0], bp_state);
    for(int i = 0; i < org_state->end_states().size(); i++) {
        org_state->end_states()[i]->get_transformed_state(bp_state, bp_state_final);
        org_state->update_end_state(i, bp_state_final);
    }
    
    math::Points t_beads (org_state->beads().size());
    math::dot_vectors(bp_state.r_T(), org_state->beads(), t_beads);
    for(int i = 0; i < t_beads.size();  i++) { t_beads[i] += bp_state.d(); }
    org_state->beads(t_beads);
    
}

inline
void
get_aligned_motif_state(
    BasepairStateOP const & ref_bp_state,
    MotifStateOP & cur_state,
    MotifStateOP const & org_state) {
    
    BasepairState bp_state, bp_state_final;
    
    ref_bp_state->get_transforming_r_and_t(*org_state->end_states()[0], bp_state);
    for(int i = 0; i < org_state->end_states().size(); i++) {
        org_state->end_states()[i]->get_transformed_state(bp_state, bp_state_final);
        cur_state->update_end_state(i, bp_state_final);
    }
    
    math::Points t_beads (org_state->beads().size());
    math::dot_vectors(bp_state.r_T(), org_state->beads(), t_beads);
    for(int i = 0; i < t_beads.size();  i++) { t_beads[i] += bp_state.d(); }
    cur_state->beads(t_beads);
    
}

#endif /* defined(__RNAMake__motif_state__) */


























