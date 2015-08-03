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
#include "secondary_structure/basepair.h"

namespace sstruct {
 
class Motif {
public:
    
    Motif():
    chains_(ChainOPs()),
    end_ids_(Strings()),
    ends_(BasepairOPs()),
    basepairs_(BasepairOPs()),
    name_(String()),
    type_("UNKNOWN")
    {}
    
    Motif(
        String const & type,
        BasepairOPs const & ends,
        ChainOPs const & chains):
    chains_(chains),
    end_ids_(Strings()),
    ends_(ends),
    basepairs_(BasepairOPs()),
    name_(String()),
    type_(type)
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
    
    inline
    BasepairOP const &
    get_bp(
        ResidueOP const & r1,
        ResidueOP const & r2) {
        
        for(auto const & bp : basepairs_) {
            if     (r1 == bp->res1() && r2 == bp->res2()) { return bp; }
            else if(r2 == bp->res1() && r1 == bp->res2()) { return bp; }
        }
        
        throw std::runtime_error("cannot find basepair in sstruct::motifs");
        
    }
    
    inline
    BasepairOP
    get_bp(
        ResidueOP const & r) {
        
        for(auto const & bp : basepairs_) {
            if(r == bp->res1() || r == bp->res2()) { return bp; }
        }
        return nullptr;
    }
    
 
public: //getters
    
    inline
    ChainOPs const &
    chains() { return chains_; }
    
    inline
    BasepairOPs const &
    basepairs() { return basepairs_; }
    
    inline
    BasepairOPs const &
    ends() { return ends_; }
    
    inline
    String const &
    type() { return type_; }
    
    inline
    String const &
    name() { return name_; }
    
    inline
    Strings const &
    end_ids() { return end_ids_; }
    
public: //setters
    
    inline
    void
    basepairs(BasepairOPs const & nbasepairs) { basepairs_ = nbasepairs; }
    
    inline
    void
    ends(BasepairOPs const & nends) { ends_ = nends; }
    
    inline
    void
    name(String const & name) { name_ = name; }
    
    inline
    void
    end_ids(Strings const & end_ids) { end_ids_ = end_ids; }
    
    
protected:
    ChainOPs chains_;
    BasepairOPs basepairs_, ends_;
    Strings end_ids_;
    String type_, name_;
    
};
    
typedef std::shared_ptr<Motif> MotifOP;
typedef std::vector<MotifOP>   MotifOPs;
    
}

#endif /* defined(__RNAMake__motif__) */
