//
//  chain.h
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__ss_chain__
#define __RNAMake__ss_chain__

#include <stdio.h>

#include "secondary_structure/residue.h"

namespace sstruct {

class Chain {
public:
    Chain() {}
    
    Chain(
        ResidueOPs residues):
    residues_(residues)
    {}
    
public:
    
    inline
    Chain
    copy() {
        ResidueOPs res;
        for (auto const & r : residues_) {
            auto r_copy = std::make_shared<Residue>(r->copy());
            res.push_back(r_copy);
        }
        return Chain(res);
    }
    
    inline
    ResidueOP const &
    first() {
        try {
            return residues_.at(0);
        }
        catch(std::out_of_range e) {
            throw SecondaryStructureException("called first() on a chain without any residues in it");
        }
        catch(...) {
            throw std::runtime_error("unexpected error in chain.first()");
        }
    }
    
    inline
    ResidueOP const &
    last() {
        try {
            return residues_.back();
        }
        catch(std::out_of_range e) {
            throw SecondaryStructureException("called last() on a chain without any residues in it");
        }
        catch(...) {
            throw std::runtime_error("unexpected error in chain.last()");
        }
    }
    
    inline
    String
    sequence() {
        String seq = "";
        for(auto const & r : residues_) { seq += r->name(); }
        return seq;
    }
    
    inline
    String
    dot_bracket() {
        String db = "";
        for(auto const & r : residues_) { db += r->dot_bracket(); }
        return db;
    }
    
    inline
    String
    to_str() {
        String s;
        for(auto const & r : residues_) { s += r->to_str() + ";"; }
        return s;
    }
    
    inline
    int
    length() {
        return (int)residues_.size();
    }
    
public: //getters
    
    inline
    ResidueOPs const &
    residues() {
        return residues_;
    }
    
private:
    ResidueOPs residues_;
    
    
};
    
Chain
str_to_chain(String const &);

    
typedef std::shared_ptr<Chain> ChainOP;
typedef std::vector<ChainOP> ChainOPs;
    
    
} //sstruct


#endif /* defined(__RNAMake__chain__) */
