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

namespace secondary_structure {

class Chain {
public:
    Chain() {}
    
    inline
    Chain(
        ResidueOPs const & residues):
    residues_(residues)
    {}
    
    inline
    Chain(
        Chain const & c) {
        residues_ = ResidueOPs(c.residues_.size());
        int i = 0;
        for(auto const & r : c.residues_) {
            residues_[i] = std::make_shared<Residue>(*r);
            i++;
        }
    }
    
    Chain(
        String const & s) {
        residues_ = ResidueOPs();
        Strings spl = base::split_str_by_delimiter(s, ";");
        for(auto const & r_str : spl) {
            if(r_str.length() < 3) { continue; }
            auto r = std::make_shared<Residue>(r_str);
            residues_.push_back(r);
        }
    }
    
public:
    
    typedef typename ResidueOPs::iterator iterator;
    typedef typename ResidueOPs::const_iterator const_iterator;
    
    iterator begin() { return residues_.begin(); }
    iterator end()   { return residues_.end(); }
    
    const_iterator begin() const { return residues_.begin(); }
    const_iterator end()   const { return residues_.end(); }
    
public:
    

    inline
    ResidueOP const &
    first() {
        if(residues_.empty()) {
            throw std::runtime_error("called first() on a chain without any residues in it");
        }

        return residues_[0];
    }
    
    inline
    ResidueOP const &
    last() {
        if(residues_.empty()) {
            throw Exception("called last() on a chain without any residues in it");
        }
     
        return residues_.back();

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

typedef std::shared_ptr<Chain> ChainOP;
typedef std::vector<ChainOP> ChainOPs;
    
    
} //secondary_structure


#endif /* defined(__RNAMake__chain__) */
