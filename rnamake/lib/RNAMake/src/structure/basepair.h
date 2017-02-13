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
#include "primitives/basepair.h"
#include "math/xyz_matrix.h"
#include "motif_state/basepair.h"
#include "structure/residue.h"


class Basepair : public primitives::Basepair {
public:
    Basepair(
            UuidOP const &,
            UuidOP const &,
            Matrix const &,
            Point const &,
            Points const &,
            StringOP const &,
            StringOP const &,
            primitives::Basepair::BasepairType const &,
            UuidOP const &);

    Basepair(
            String const &,
            UuidOP const &,
            UuidOP const &);

    Basepair(Basepair const &);

    Basepair(
            Basepair const & bp,
            UuidOP const & res1_uuid,
            UuidOP const & res2_uuid,
            UuidOP const & bp_uuid = nullptr);


    ~Basepair() {
        res1_uuid_.reset();
        res2_uuid_.reset();
    }

public:

    void
    move(Point const & p) {
        sugars_[0] += p;
        sugars_[1] += p;
        d_ += p;
    }

    inline
    void
    transform(Transform const & t) {
        auto r_T = t.rotation().transpose();
        auto new_r = Matrix();
        dot(r_, r_T, new_r);
        r_ = new_r;
        sugars_[0] = dot_vector(r_T, sugars_[0]) + t.translation();
        sugars_[1] = dot_vector(r_T, sugars_[1]) + t.translation();
        d_ = dot_vector(r_T, d_) + t.translation();
    }

    inline
    void
    fast_transform(
            Matrix const & r,
            Vector const & t) {

    }

    inline
    float
    diff(Basepair const & bp) {
        return 0.0f;
    }

private:

    inline
    float
    _rot_diff(Basepair const & bp) {
        return 0.0f;
    }

public:

    void
    flip_res();

    String
    to_str() const;

    state::BasepairOP
    get_state();


public: // getters

    inline
    Matrix const &
    r() const { return r_; }

    inline
    Point const &
    d() const { return d_; }

    inline
    Points const &
    sugars() const { return sugars_; }

    inline
    Point const &
    res1_sugar() const { return sugars_[0]; }

    inline
    Point const &
    res2_sugar() const { return sugars_[1]; }

    inline
    StringOP const &
    x3dna_bp_type() const { return x3dna_bp_type_; }

    inline
    UuidOP const &
    uuid() const { return uuid_; }

    inline
    UuidOP const &
    res1_uuid() const { return res1_uuid_; }

    inline
    UuidOP const &
    res2_uuid() const { return res2_uuid_; }

    inline
    StringOP const &
    name() const { return name_; }

    inline
    primitives::Basepair::BasepairType
    bp_type() const { return bp_type_; }

    /*inline
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
    
    */
private:
    Matrix r_;
    Point d_;
    Points sugars_;
    UuidOP res1_uuid_, res2_uuid_;
    StringOP name_, x3dna_bp_type_;
    primitives::Basepair::BasepairType bp_type_;

};

typedef std::shared_ptr<Basepair> BasepairOP;
typedef std::vector<BasepairOP> BasepairOPs;

bool
are_basepairs_equal(
        BasepairOP const & bp1,
        BasepairOP const & bp2,
        int check_uuids = 1);

Point
_calc_center(ResidueOPs const &);



#endif /* defined(__RNAMake__basepair__) */
