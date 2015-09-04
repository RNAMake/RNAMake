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
        String const &);
    
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
    
private:
    String id_;
    int block_end_add_;
    MotifStateEnsembleMemberOPs members_;
    RandomNumberGenerator rng_;
    
};

typedef std::shared_ptr<MotifStateEnsemble> MotifStateEnsembleOP;
typedef std::vector<MotifStateEnsembleOP>   MotifStateEnsembleOPs;


#endif /* defined(__RNAMake__motif_state_ensemble__) */






















