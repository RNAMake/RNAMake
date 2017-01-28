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

//RNAMake Headers
#include "base/types.h"
#include "util/uuid.h"
#include "math/xyz_matrix.h"
#include "structure/atom.h"
#include "structure/residue.h"
#include "structure/basepair_state.h"

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
        /*for( auto const & a : res1->atoms() ) {
            if(a != nullptr) { atoms_.push_back(a); }
        }
        for( auto const & a : res2->atoms() ) {
            if(a != nullptr) { atoms_.push_back(a); }
        }*/
        
        Point d = calc_center(atoms_);
        Points sugars(2);
        sugars[0] = res1_->get_atom("C1'")->coords();
        sugars[1] = res2_->get_atom("C1'")->coords();
        bp_state_ = std::make_shared<BasepairState>(d, r, sugars);
        
    }
    
    ~Basepair() {
        res1_.reset();
        res2_.reset();
    }
    
    Basepair
    copy();
    
    inline
    bool
    operator == (Basepair const & bp ) const { return uuid_ == bp.uuid_; }

public: // getters
    
    inline
    BasepairStateOP const &
    state() {
        bp_state_->d(calc_center(atoms_));
        bp_state_->sugars(Points{ res1_->get_atom("C1'")->coords(), res2_->get_atom("C1'")->coords() });
        return bp_state_;
    }
    
    inline
    float
    diff(Basepair & bp) {
        auto diff = d().distance(bp.d());
        diff += _rot_diff(bp)*2;
        return diff;
        
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
        if     ( res->uuid() == res1_->uuid()) {  return res2_;  }
        else if( res->uuid() == res2_->uuid()) {  return res1_;  }
        else { throw "called partner with resiude not in this basepair"; }
    }
    
    inline
    String const
    name() const {
        auto res1_name = res1_->chain_id()+std::to_string(res1_->num())+res1_->i_code();
        auto res2_name = res2_->chain_id()+std::to_string(res2_->num())+res2_->i_code();
        
        if(res1_->chain_id() < res2_->chain_id()) { return res1_name+"-"+res2_name; }
        if(res1_->chain_id() > res2_->chain_id()) { return res2_name+"-"+res1_name; }
        
        if(res1_->num() < res2_->num()) { return res1_name+"-"+res2_name; }
        else                            { return res2_name+"-"+res1_name; }
        
        return res1_name+"-"+res2_name;
    }
    
    inline
    void
    flip() { bp_state_->flip(); }
    
    inline
    Matrix const &
    r() const { return bp_state_->r(); }
    
    inline
    Point const
    d()  const { return calc_center(atoms_); }
    
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
    
    inline
    AtomOPs const &
    atoms() const { return atoms_; }

public: // setters
    
    inline
    void
    r(Matrix const & nr) { bp_state_->r(nr); }
    
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
    
    
    inline
    void
    bp_type(String const & nbp_type) { bp_type_ = nbp_type; }
    
public:
    
    String const
    to_str() const;
    
    String const
    to_pdb_str() const;
    
    void
    to_pdb(String const) const;
    
private:
    inline
    float
    _rot_diff(
        Basepair & bp) {
        auto r_diff = r().difference(bp.r());
        bp.flip();
        float r_diff_2 = r().difference(bp.r());
        bp.flip();
        
        if( r_diff > r_diff_2) { r_diff = r_diff_2;}
        
        return r_diff;
    }
    
    
private:
    ResidueOP res1_, res2_;
    AtomOPs atoms_;
    BasepairStateOP bp_state_;
    String bp_type_;
    int flipped_;
    Uuid uuid_;
    
};

typedef std::shared_ptr<Basepair> BasepairOP;
typedef std::vector<BasepairOP> BasepairOPs;

bool
wc_bp(BasepairOP const &);

bool
gu_bp(BasepairOP const &);


#endif /* defined(__RNAMake__basepair__) */
