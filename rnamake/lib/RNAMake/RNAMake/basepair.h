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
    Basepair()
    {}
    
    Basepair(
        ResidueOP & res1,
        ResidueOP & res2,
        Matrix const & r,
        String const & bp_type):
        bp_type_(bp_type),
        uuid_( Uuid() ),
        res1_( res1 ),
        res2_( res2 ),
        flipped_ ( 0 )
    {
    
        atoms_ = AtomOPs();
        for( auto const & a : res1->atoms() ) {
            if(a != NULL) { atoms_.push_back(a); }
        }
        for( auto const & a : res2->atoms() ) {
            if(a != NULL) { atoms_.push_back(a); }
        }
        
        Point d = center(atoms_);
        Points sugars(2);
        sugars[0] = res1_->get_atom("C1'")->coords();
        sugars[1] = res2_->get_atom("C1'")->coords();
        bp_state_ = BasepairState(d, r, sugars);
        
    }
    
    ~Basepair() {}
    
    Basepair
    copy();
    
    inline
    bool
    operator == (Basepair const & bp ) const { return uuid_ == bp.uuid_; }

public: // getters
    
    inline
    BasepairState &
    state() {
        bp_state_.d(center(atoms_));
        return bp_state_;
    }
    
    inline
    ResidueOPs const
    residues() const {
        ResidueOPs res(2);
        res[0] = res1_;
        res[1] = res2_;
        return res;
    }
    
    inline
    ResidueOP const &
    partner(ResidueOP const & res) {
        if     ( res == res1_) {  return res2_;  }
        else if( res == res2_) {  return res1_;  }
        else { throw "called partner with resiude not in this basepair"; }
    }
    
    inline
    String const
    name() const {
        std::stringstream ss;
        ss << res1_->chain_id() << res1_->num() << res1_->i_code();
        ss << "-";
        ss << res2_->chain_id() << res2_->num() << res2_->i_code();
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
    Point const
    d()  const { return center(atoms_); }
    
    inline
    Uuid const &
    uuid() const { return uuid_; }
    
    inline
    ResidueOP
    res1() const { return res1_; }
    
    inline
    ResidueOP 
    res2() const { return res2_; }
    
    inline
    String const &
    bp_type() const { return bp_type_; }
    
    inline
    int const
    flipped() const { return flipped_; }

public: // setters
    
    inline
    void
    r(Matrix const & nr) { bp_state_.r(nr); }
    
    inline
    void
    uuid(Uuid const & nuuid) { uuid_ = nuuid; }
    
    inline
    void
    flipped(int const & nflipped) { flipped_ = nflipped; }
    
    inline
    void
    res1(ResidueOP const & nres1) { res1_ = nres1;}
    
    inline
    void
    res2(ResidueOP const & nres2) { res2_ = nres2;}
    
public:
    
    String const
    to_str() const;
    
    String const
    to_pdb_str() const;
    
    void
    to_pdb(String const) const;
    
    
private:
    ResidueOP res1_, res2_;
    AtomOPs atoms_;
    BasepairState bp_state_;
    String bp_type_;
    int flipped_;
    Uuid uuid_;
    
};

typedef std::vector<Basepair> Basepairs;
typedef std::shared_ptr<Basepair> BasepairOP;
typedef std::vector<BasepairOP> BasepairOPs;

bool
wc_bp(BasepairOP const &);

bool
gu_bp(BasepairOP const &);


#endif /* defined(__RNAMake__basepair__) */
