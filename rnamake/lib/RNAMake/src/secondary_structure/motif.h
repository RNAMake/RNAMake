//
//  motif.h
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__ss_motif__
#define __RNAMake__ss_motif__

#include <stdio.h>

#include "secondary_structure/chain.h"

namespace sstruct {
 
class Motif {
public:
    
    Motif():
    chains_(ChainOPs()),
    end_ids_(Strings())
    {}
    
    ~Motif() {}
    
    
public:
    
    inline
    ResidueOPs
    residues() {
        ResidueOPs res;
        for(auto const & c : chains_) {
            for(auto const & r : c->residues()) {
                res.push_back(r);
            }
        }
        return res;
    }
    
    inline
    String
    sequence() {
        String seq;
        int i = 0;
        for (auto const & c : chains_) {
            seq += c->sequence();
            if(i+1 != chains_.size()) { seq += "&"; }
            i++;
        }
        return seq;
    }
    
    inline
    String
    dot_bracket() {
        String db;
        int i = 0;
        for (auto const & c : chains_) {
            db += c->dot_bracket();
            if(i+1 != chains_.size()) { db += "&"; }
            i++;
        }
        return db;
    }
    
    inline
    ResidueOP
    get_residue(
        int num,
        String const & chain_id,
        String i_code = "") {
        
        for(auto & c : chains_) {
            for (auto & r : c->residues() ){
                if (num == r->num() && chain_id == r->chain_id() && i_code == r->i_code()) {
                    return r;
                }
            }
        }
        return nullptr;
        
    }

    inline
    ResidueOP
    get_residue(
        Uuid const & uuid) {
        
        for(auto & c : chains_) {
            for (auto & r : c->residues() ){
                if (r->uuid() == uuid) {
                    return r;
                }
            }
        }
        return nullptr;

    }
    
public: //getters
    
    inline
    ChainOPs
    chains() { return chains_; }
    
    
protected:
    ChainOPs chains_;
    Strings end_ids_;
    
    
};
    
}

#endif /* defined(__RNAMake__motif__) */
