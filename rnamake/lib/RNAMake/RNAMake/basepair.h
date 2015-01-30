//
//  basepair.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__basepair__
#define __RNAMake__basepair__

#include <stdio.h>
#include "basepair_state.h"
#include "residue.h"
#include "atom.h"
#include "uuid.h"
#include "xyzMatrix.h"
#include "types.h"
#include "util.h"

class Basepair {
public:
    Basepair() {}
    
    Basepair(
        Residue const & res1,
        Residue const & res2,
        Matrix const & r,
        String const & bp_type):
        res1_(res1),
        res2_(res2),
        bp_type_(bp_type),
        uuid_( Uuid() ),
        flipped_ ( 0 )
    {
        
        atoms_ = AtomOPs();
        for( auto const & a : res1.atoms() ) {
            if(a != NULL) { atoms_.push_back(a); }
        }
        for( auto const & a : res2.atoms() ) {
            if(a != NULL) { atoms_.push_back(a); }
        }
        
        Point d = center(atoms_);
        Points sugars(2);
        sugars[0] = res1_.get_atom("C1'")->coords();
        sugars[1] = res2_.get_atom("C1'")->coords();
        bp_state_ = BasepairState(d, r, sugars);
        
    }
    
    ~Basepair() {}
    
    Basepair
    copy();

public: // getters
    
    inline
    BasepairState const &
    state() {
        bp_state_.d(center(atoms_));
        return bp_state_;
    }

    inline
    Residues const
    residues() {
        Residues res(2);
        res[0] = res1_;
        res[1] = res2_;
        return res;
    }
    
    inline
    Residue const &
    partner(Residue const & res) {
        if     ( res == res1_) {  return res2_;  }
        else if( res == res2_) {  return res1_;  }
        else { throw "called partner with resiude not in this basepair"; }
    }
    
    inline
    String const
    name() {
        std::stringstream ss;
        ss << res1_.short_name() << res1_.num() << res1_.i_code();
        ss << res2_.short_name() << res2_.num() << res2_.i_code();
        return ss.str();
    }
    
    inline
    void
    flip(int flip_d=-1) {
        if( flip_d == -1) {
            if(flipped_ == 0) {  flip_d = 1; }
            else              {  flip_d = 0; }
            bp_state_.flip();
            flipped_ = flip_d;
            return;
        }
        
        if( flip_d == flipped_) { return; }
        else {
            bp_state_.flip();
            flipped_ = flip_d;
        }
    }
    
    inline
    Matrix const &
    r() const { return bp_state_.r(); }
    
    inline
    Point const &
    d()  { return state().d(); }
    
public:
    
    String const
    to_str();
    
    String const
    to_pdb_str();
    
    void
    to_pdb(String const);
    
    
private:
    Residue res1_, res2_;
    AtomOPs atoms_;
    BasepairState bp_state_;
    String bp_type_;
    int flipped_;
    Uuid uuid_;
    
};

typedef std::vector<Basepair> Basepairs;

#endif /* defined(__RNAMake__basepair__) */
