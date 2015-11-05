//
//  motif_state_ensemble.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_ensemble__
#define __RNAMake__motif_state_ensemble__

#include <stdio.h>
#include <algorithm>

#include "util/random_number_generator.h"
#include "motif/motif_state.h"

struct MotifStateEnsembleMember {
    inline
    MotifStateEnsembleMember(
        MotifStateOP const & nmotif_state,
        float const & nenergy):
    motif_state(nmotif_state),
    energy(nenergy)
    {}
    
    inline
    MotifStateEnsembleMember
    copy() {
        return MotifStateEnsembleMember(std::make_shared<MotifState>(motif_state->copy()),
                                        energy);
    }
    
    inline
    String
    to_str() {
        return motif_state->to_str() + "#" + std::to_string(energy);
    }
    
    MotifStateOP motif_state;
    float energy;
};

typedef std::shared_ptr<MotifStateEnsembleMember> MotifStateEnsembleMemberOP;
typedef std::vector<MotifStateEnsembleMemberOP>   MotifStateEnsembleMemberOPs;

//to allow sorting of members
struct MotifStateEnsembleMember_LessThanKey {
    inline
    bool
    operator() (
        MotifStateEnsembleMemberOP const & mem1,
        MotifStateEnsembleMemberOP const & mem2) {
        return mem1->energy < mem2->energy;
    }
};


class MotifStateEnsemble {
public:
    MotifStateEnsemble():
    id_(""),
    block_end_add_(0),
    members_(MotifStateEnsembleMemberOPs()),
    rng_(RandomNumberGenerator())
    {}
    
    MotifStateEnsemble(
        MotifStateOP const & ms):
    id_(ms->end_ids()[0]),
    block_end_add_(ms->block_end_add()),
    members_(MotifStateEnsembleMemberOPs()),
    rng_(RandomNumberGenerator()) {
        
        members_.push_back(std::make_shared<MotifStateEnsembleMember>(ms, 1));
    }
    
    MotifStateEnsemble(
        String const &);
    
    ~MotifStateEnsemble() {}
    
public:
    
    void
    setup(
        String const &,
        MotifStateOPs const &,
        Floats const &);
    
    MotifStateEnsemble
    copy();
    
    String
    to_str();
    
    inline
    MotifStateEnsembleMemberOP const &
    get_random_member() {
        return members_[rng_.randrange((int)members_.size()-1) ];
    }
    
    inline
    int
    num_end_states() { return (int)members_[0]->motif_state->end_states().size(); }
    
    inline
    MotifStateOP const &
    most_populated() { return members_[0]->motif_state; }
    
    inline
    MotifStateEnsembleMemberOP const &
    get_member(
        int i) {
        return members_[i];
    }
    
    inline
    int
    member_index(
        MotifStateEnsembleMemberOP const & mem) {
        return std::find(members_.begin(), members_.end(), mem) - members_.begin();
        
    }
    
public:
    
    inline
    int
    block_end_add() { return block_end_add_; }
    
    inline
    MotifStateEnsembleMemberOPs members() { return members_; }
    
    inline
    String
    id() { return id_; }
    
private:
    String id_;
    int block_end_add_;
    MotifStateEnsembleMemberOPs members_;
    RandomNumberGenerator rng_;
    
};

typedef std::shared_ptr<MotifStateEnsemble> MotifStateEnsembleOP;
typedef std::vector<MotifStateEnsembleOP>   MotifStateEnsembleOPs;


#endif /* defined(__RNAMake__motif_state_ensemble__) */






















