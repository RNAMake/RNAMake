//
//  structure.h
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__structure__
#define __RNAMake__structure__

#include <stdio.h>

//RNAMake Headers
#include "base/types.h"
#include "secondary_structure/chain.h"

namespace sstruct {
    
class Structure {
public:
    inline
    Structure(
        ChainOPs const & chains):
    chains_(chains)
    {}
    
    Structure(
        String const & sequence,
        String const & dot_bracket):
    chains_(ChainOPs()){
        _setup_chains(sequence, dot_bracket);
    }
    
    Structure(
        Structure const & structure) {
    
        chains_ = ChainOPs(structure.chains_.size());
        int i = 0;
        for(auto const & c : structure.chains_) {
            chains_[i] = std::make_shared<Chain>(*c);
            i++;
        }
    }
    
public:
    
    inline
    ResidueOPs
    residues() {
        ResidueOPs res;
        for (auto const c : chains_) {
            std::copy(c->residues().begin(),
                      c->residues().end(),
                      std::inserter(res, res.end()));
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
    
    ResidueOP const
    get_residue(
        int const & ,
        String const & ,
        String const & );
    
    ResidueOP const
    get_residue(
        Uuid const &);
    
    String
    to_str() {
        String s = "";
        for(auto const & c : chains_) {
            s += c->to_str() + "|";
        }
        return s;
    }
    
public:
    inline
    ChainOPs const &
    chains() { return chains_; }


private:
    
    void
    _setup_chains(
        String const &,
        String const &);
    
private:
    ChainOPs chains_;
};
   

Structure
str_to_structure(
    String const &);
    

typedef std::shared_ptr<Structure> StructureOP;

}




#endif /* defined(__RNAMake__structure__) */