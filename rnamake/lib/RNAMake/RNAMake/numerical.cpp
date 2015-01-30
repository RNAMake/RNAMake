//
//  Numeric_Test.cpp
//  REDESIGNC
//
//  Created by Joseph Yesselman on 10/4/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#include "numerical.h"
#include "xyzVector.h"
#include "xyzMatrix.h"
#include "basepair_state.h"

int
are_floats_equal(
	float const a,
	float const b) {
	if (fabs(a - b) < 0.001) { return 1; }
	else					 { return 0; }
}

int
are_xyzVector_equal(
	Vector const & vec,
	Vector const & correct_vec) {
    
	if (are_floats_equal(vec.x(),correct_vec.x()) &&
		are_floats_equal(vec.y(),correct_vec.y()) &&
		are_floats_equal(vec.z(),correct_vec.z())) {
		return 1;
	}
	else {
		return 0;
	}
}

int
are_xyzVectors_equal(
	Vectors const & v,
	Vectors const & vc) {
	
	if(v.size() != vc.size()) { return 0; }
	
	for(int i = 0; i < v.size(); i++) {
		if(! are_xyzVector_equal(v[i], vc[i])) { return 0; }
	}
	
	return 1;
	
}


int
are_xyzMatrix_equal(
	Matrix const & m,
	Matrix const & mc) {
	
	if(!are_floats_equal(m.xx(),mc.xx()) ||
	   !are_floats_equal(m.xy(),mc.xy()) ||
	   !are_floats_equal(m.xz(),mc.xz()) ||
	   !are_floats_equal(m.yx(),mc.yx()) ||
	   !are_floats_equal(m.yz(),mc.yz()) ||
	   !are_floats_equal(m.zx(),mc.zx()) ||
	   !are_floats_equal(m.zy(),mc.zy()) ||
	   !are_floats_equal(m.zz(),mc.zz())) {
		return 0;
	}
	
	return 1;
	
}


int
are_BasepairStates_equal(
	BasepairState const & a,
	BasepairState const & b) {
	
	if(!are_xyzVector_equal(a.d(),b.d()) ||
	   !are_xyzMatrix_equal(a.r(),b.r())) {
		return 0;
	}
	
	for(int i = 0; i < 2; i++) {
		if(!are_xyzVector_equal(a.sugars()[i], b.sugars()[i])) { return 0; }
	}
	
	return 1;
}
