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

/*
 * Exception for motif ensemble
 */
class MotifEnsembleException : public std::runtime_error {
public:
    /**
     * Standard constructor for MotifEnsembleException
     * @param   message   Error message for motif ensembles
     */
    MotifEnsembleException(String const & message):
    std::runtime_error(message)
    {}
};


class MotifEnsemble {
private:
    struct MotifEnsembleMember {
        MotifEnsembleMember(
            MotifOP const & nmotif,
            float const & nenergy):
        motif(nmotif),
        energy(nenergy)
        {}
        
        String
        to_str() { return motif->to_str() + "#" + std::to_string(energy); }
        
        MotifOP motif;
        float energy;
    };
    
    
    typedef std::shared_ptr<MotifEnsembleMember> MotifEnsembleMemberOP;
    typedef std::vector<MotifEnsembleMemberOP>   MotifEnsembleMemberOPs;


public:
    MotifEnsemble():
    id_(""),
    block_end_add_(0),
    members_(MotifEnsembleMemberOPs())
    {}
    
    MotifEnsemble(
        MotifOPs const & motifs,
        Floats const & energies):
    members_(MotifEnsembleMemberOPs()),
    block_end_add_(0) {
        
        if(motifs.size() != energies.size()) {
            throw MotifEnsembleException(
                "needs to be an equal number of motifs and energies");
        }
        
        if(motifs.size() > 0) {
            block_end_add_ = motifs[0]->block_end_add();
            
            for(auto const & m : motifs) {
                if(block_end_add_ != m->block_end_add()) {
                    throw MotifEnsembleException(
                        "all motifs included in the ensemble must have the same block_end_add_");
                }
            }
            
        }
        
        int i = 0;
        for(auto const & m : motifs) {
            members_.push_back(std::make_shared<MotifEnsembleMember>(m, energies[i]));
        }
    }
    
    MotifEnsemble(
        MotifEnsemble const & me):
    id_(me.id_),
    block_end_add_(me.block_end_add_),
    members_(MotifEnsembleMemberOPs()) {
        for(auto const & mem : me.members_) {
            auto m = std::make_shared<Motif>(*mem->motif);
            members_.push_back(std::make_shared<MotifEnsembleMember>(m, mem->energy));
        }
    }
    
    MotifEnsemble(
        String const & s,
        ResidueTypeSet const & rts) {
        
        auto spl = split_str_by_delimiter(s, "{");
        id_ = spl[0];
        block_end_add_ = std::stoi(spl[1]);
        members_ = MotifEnsembleMemberOPs();
        for(int i = 2; i < spl.size()-1; i++) {
            auto m_spl = split_str_by_delimiter(spl[i], "#");
            auto m = std::make_shared<Motif>(m_spl[0], rts);
            auto energy = std::stof(m_spl[1]);
            members_.push_back(std::make_shared<MotifEnsembleMember>(m, energy));
        }
        
    }
    
    
    ~MotifEnsemble() {}
    
public:
    
    size_t
    size() { return members_.size(); }
    
    String
    to_str() {
        auto s = String(id_ + "{" + std::to_string(block_end_add_) + "{");
        for(auto const & mem : members_ ) {
            s += mem->to_str() + "{";
        }
        
        return s;
    }
    
    MotifStateEnsembleOP
    get_state();
    
    MotifEnsembleMemberOPs const &
    members() { return members_; }
    
    
private:
    String id_;
    int block_end_add_;
    MotifEnsembleMemberOPs members_;
    
};

typedef std::shared_ptr<MotifEnsemble> MotifEnsembleOP;



#endif /* defined(__RNAMake__motif_ensemble__) */
