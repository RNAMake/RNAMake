//
//  BasepairState.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 9/28/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __REDESIGNC__BasepairState__
#define __REDESIGNC__BasepairState__

#include <iostream>
#include <vector>

//RNAMake Headers
#include "math/xyz_vector.h"
#include "math/xyz_matrix.h"
#include "math/transform.h"
#include "math/numerical.h"
#include "structure/basepair_state.fwd.h"

class BasepairState {
	
public:
	
	inline
	BasepairState():
	d_( Point(0.0) ),
	r_( Matrix(0.0) ),
	r_T_( Matrix(0.0) ),
	sugars_ ( Points(2) ),
	diff_ ( Vector(0.0) ),
	diff_sugars_( Vectors(2) )
 	{}
	
	inline
	BasepairState(
		Point const & d,
		Matrix const & r,
		Points const & sugars):
	d_( d ),
	r_( r ),
	r_T_( Matrix(0) ),
	sugars_( sugars ),
	diff_ ( Vector(0.0) ),
	diff_sugars_( Vectors(2) )
	{ transpose(r_, r_T_); }
	
    inline
    BasepairState(BasepairState const & b):
    d_( b.d_ ),
    r_( b.r_ ),
    r_T_( Matrix(0) ),
    sugars_ (b.sugars_),
    diff_ ( Vector(0.0) ),
    diff_sugars_( Vectors(2) )
    { transpose(r_, r_T_); }

    
	~BasepairState()
	{}
    
    BasepairState
    copy() const {
        //Point d = d_;
        //Matrix r = r_;
        //Points sugars = sugars_;
        return BasepairState(d_, r_, sugars_);
    }

	
public:
	inline
	void
	calculate_r_T() { transpose(r_, r_T_); }
	
	inline
	void
	get_transforming_r_and_t (
		BasepairState const & o_state, //state with desired rotation and translation
		BasepairState & r_state) {

        //TODO figure out why I have to do this each time???
        calculate_r_T();
        
		//calculate transforming rotation matrix and store it in r_state (result state)
		dot(r_T_,o_state.r_,r_state.r_);
		r_state.calculate_r_T();
		
		//rotate sugars to new position and store them in r_state
		dot_vectors(r_state.r_T_,o_state.sugars_,r_state.sugars_);
		
		diff_ = -o_state.d() + d_;
		int i;
		
		for(i = 0; i < 2; i+=1) {
			r_state.sugars_[i] += diff_;
		}
		
		for(i = 0; i < 2; i+=1) {
			diff_sugars_[i] = sugars_[i] - r_state.sugars_[i];
		}
		
		diff_ = (diff_sugars_[0] + diff_sugars_[1]) / 2.0f;
		r_state.d_ = -o_state.d() + diff_ + d_;
		
	}
	
    inline
	void
	get_transformed_state(
		BasepairState const & o_state,
		BasepairState & r_state) {
		
		dot        (r_, o_state.r_T_, r_state.r_);
		dot_vector (o_state.r_T_, d_, r_state.d_);
		dot_vectors(o_state.r_T_, sugars_, r_state.sugars_);

		int i;
		for(i = 0; i < 2; i++) {
			r_state.sugars_[i] += o_state.d_;
		}
		r_state.d_ += o_state.d_;
		
	}
	
    inline
    void
    flip() {
        r_ = transform_1(r_);
        calculate_r_T();
    }
    
    inline
    float
    diff(
        BasepairStateOP const & state) {
        float diff = d_.distance(state->d());
        diff += _rot_diff(state);
        return diff;
    }
    
    inline
    float
    _rot_diff(
        BasepairStateOP const & state) {
        float r_diff = r_.difference(state->r());
        state->flip();
        float r_diff_2 = r_.difference(state->r());
        state->flip();
        
        if( r_diff > r_diff_2) { r_diff = r_diff_2;}
        
        return r_diff;
    }
    
public:
    
    String const
    to_str() const{
        String s = vector_to_str(d_) + ";" + matrix_to_str(r_) + ";" + vector_to_str(sugars_[0]) + vector_to_str(sugars_[1]);
        return s;
    }
	
	
public: //getters
	
	inline
	const
	Point &
	d() {
		return d_;
	}
	
	inline
	const
	Point &
	d() const {
		return d_;
	}
	
	inline
	const
	Matrix &
	r() {
		return r_;
	}
	
	inline
	const
	Matrix &
	r() const {
		return r_;
	}
	
	inline
	const
	Matrix &
	r_T() {
		return r_T_;
	}
	
	inline
	const
	Matrix &
	r_T() const {
		return r_T_;
	}
	
	inline
	const
	Points &
	sugars() {
		return sugars_;
	}
	
	inline
	const
	Points &
	sugars() const {
		return sugars_;
	}


public: //setters
	
	inline
	void
	d(Point const & newd) { d_ = newd; }
	
	inline
	void
	r( Matrix const & newr) { r_ = newr; }
	
	inline
	void
	sugars( Points const & newsug) { sugars_ = newsug; }
	
    inline
	void
	set (BasepairState const & nstate_) {
		d_ = nstate_.d_;
		r_ = nstate_.r_;
		sugars_ = nstate_.sugars_;
		calculate_r_T();
	}
	
	
private:	
	Point d_;
	Matrix r_;
	Matrix r_T_; //holds the transpose of r_ for speed
	Points sugars_;
	
	/*
	Inclusion of these two tempory variables increases the speed of
	get_transforming_r_and_t() by 240 percent seems worth it for just 
	a few more bytes of memory since this a bulk of the calculations 
	done in the algorithm -JDY 2014.10.4
	*/
	
	Vector diff_; //stores partial products to speed up algorithm
	Vectors diff_sugars_; //stores partial products to speed up algorithm

	
};

BasepairState
str_to_basepairstate(
	String const &);

BasepairState
get_ref_bp_state();

float
get_bpstate_rotation_diff(
	BasepairState const &,
	BasepairState const &);


std::ostream&
operator <<(
	std::ostream&,
	const BasepairState&);


inline
const
float
frame_distance(
	BasepairStateOP const & current,
	BasepairStateOP const & end,
	BasepairStateOP const & endflip) {
	
	float score = current->d().distance(end->d());
	
	float r_diff      = end->r().difference(current->r());
	float r_diff_flip = endflip->r().difference(current->r());
	
	if(r_diff > r_diff_flip) {
		r_diff = r_diff_flip;
	}
	
	score += r_diff;
	
	return score;
	
}

inline
const
float
new_score_function(
	BasepairStateOP const & current,
	BasepairStateOP const & end,
	BasepairStateOP const & endflip) {
	
	float d_diff = current->d().distance(end->d());
	
	if(d_diff > 25) { return d_diff; }
	
	float r_diff      = current->r().difference(end->r());
	float r_diff_flip = current->r().difference(endflip->r());
	
	if(r_diff > r_diff_flip) {
		r_diff = r_diff_flip;
	}
	
	float scale = (log(150/d_diff) - 1);
	if (scale > 2) { scale = 2; }
	if (scale < 0) { scale = 0; }
	
	return d_diff + scale*r_diff;
}

inline
const
float
new_score_function_new(
    BasepairStateOP const & current,
    BasepairStateOP const & end,
    BasepairStateOP const & endflip) {
    
    float d_diff = (current->sugars()[0].distance(end->sugars()[1]) +
                    current->sugars()[1].distance(end->sugars()[0]))*0.50;
    
    if(d_diff > 25) { return d_diff; }
    
    float r_diff      = current->r().difference(end->r());
    float r_diff_flip = current->r().difference(endflip->r());
    
    if(r_diff > r_diff_flip) {
        r_diff = r_diff_flip;
    }
    
    float scale = (log(150/d_diff) - 1);
    if (scale > 2) { scale = 2; }
    if (scale < 0) { scale = 0; }
    
    return d_diff + scale*r_diff;
}


int
are_BasepairStates_equal(
    BasepairState const,
    BasepairState const);



#endif /* defined(__REDESIGNC__BasepairState__) */
