//
//  basepair.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__STRUCTURE_basepair__
#define __RNAMake__STRUCTURE_basepair__

#include <stdio.h>
#include <util/x3dna.h>

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
            Uuid const &,
            Uuid const &,
            Matrix const &,
            Point const &,
            Point const &,
            Point const &,
            SimpleStringOP const &,
            X3dna::X3dnaBPType const,
            primitives::Basepair::BasepairType const,
            Uuid const &);

    Basepair(
            String const &,
            Uuid const &,
            Uuid const &);

    Basepair(Basepair const &);

    Basepair(
            Basepair const & bp,
            Uuid const &,
            Uuid const &,
            Uuid const & );

    ~Basepair() {}

public:

    void
    move(Point const & p) {
        res1_sugar_ += p;
        res2_sugar_ += p;
        d_ += p;
    }

    inline
    void
    transform(Transform const & t) {
        auto r_T = t.rotation().transpose();
        auto new_r = Matrix();
        dot(r_, r_T, new_r);
        r_ = new_r;
        res1_sugar_ = dot_vector(r_T, res1_sugar_) + t.translation();
        res2_sugar_ = dot_vector(r_T, res2_sugar_) + t.translation();
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
    Points
    sugars() const { return Points{res1_sugar_, res2_sugar_}; }

    inline
    Point const &
    res1_sugar() const { return res1_sugar_; }

    inline
    Point const &
    res2_sugar() const { return res2_sugar_; }

    inline
    X3dna::X3dnaBPType
    x3dna_bp_type() const { return x3dna_bp_type_; }

    inline
    Uuid const &
    uuid() const { return uuid_; }

    inline
    SimpleString const &
    name() const { return *name_; }

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
    ResidueOP const &
    partner(ResidueOP const & res) {
        if     ( res->uuid() == res1_->uuid()) {  return res2_;  }
        else if( res->uuid() == res2_->uuid()) {  return res1_;  }
        else { throw "called partner with resiude not in this basepair"; }
    }



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
    Point res1_sugar_, res2_sugar_;
    SimpleStringOP name_;
    X3dna::X3dnaBPType x3dna_bp_type_;
    primitives::Basepair::BasepairType bp_type_;

};

typedef std::shared_ptr<Basepair> BasepairOP;

bool
are_basepairs_equal(
        BasepairOP const & bp1,
        BasepairOP const & bp2,
        int check_uuids = 1);

typedef std::vector<BasepairOP> BasepairOPs;

Point
_calc_center(ResidueOPs const &);



#endif /* defined(__RNAMake__STRUCTURE_basepair__) */
