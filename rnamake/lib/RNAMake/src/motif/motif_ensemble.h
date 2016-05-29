//
//  motif_ensemble.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_ensemble__
#define __RNAMake__motif_ensemble__

#include <stdio.h>

//RNAMake Headers
#include "motif/motif.h"
#include "motif/motif_state_ensemble.h"

struct MotifEnsembleMember {
    MotifEnsembleMember(
        MotifOP const & nmotif,
        float const & nenergy):
    motif(nmotif),
    energy(nenergy)
    {}
    
    MotifOP motif;
    float energy;
};

typedef std::shared_ptr<MotifEnsembleMember> MotifEnsembleMemberOP;
typedef std::vector<MotifEnsembleMemberOP>   MotifEnsembleMemberOPs;

class MotifEnsemble {
public:
    MotifEnsemble():
    id_(""),
    block_end_add_(0),
    members_(MotifEnsembleMemberOPs())
    {}
    
    MotifEnsemble(
        String const & s):
    id_(""),
    block_end_add_(0),
    members_(MotifEnsembleMemberOPs()) {
        
        ResidueTypeSet rts;
        auto spl = split_str_by_delimiter(s, "{");
        id_ = spl[0];
        block_end_add_ = std::stoi(spl[1]);
        for(int i = 2; i < spl.size(); i++) {
            auto spl2 = split_str_by_delimiter(spl[i], "#");
            auto m = std::make_shared<Motif>(spl2[0], rts);
            auto energy = std::stof(spl2[1]);
            auto ms = std::make_shared<MotifEnsembleMember>(m, energy);
            members_.push_back(ms);
        }
    }
    
    MotifStateEnsembleOP
    get_state() {
        auto mse = std::make_shared<MotifStateEnsemble>();
        auto motif_states = MotifStateOPs();
        auto energies = Floats();
        
        for(auto const & mem : members_) {
            motif_states.push_back(mem->motif->get_state());
            energies.push_back(mem->energy);
        }
        mse->setup(id_, motif_states, energies);
        return mse;
    }

public:
    inline
    MotifEnsembleMemberOPs const &
    members() { return members_; }

private:
    String id_;
    int block_end_add_;
    MotifEnsembleMemberOPs members_;
    
};

typedef std::shared_ptr<MotifEnsemble> MotifEnsembleOP;



#endif /* defined(__RNAMake__motif_ensemble__) */
