//
//  Numeric_Test.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 10/4/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __REDESIGNC__Numeric_Test__
#define __REDESIGNC__Numeric_Test__

#include <iostream>

//RNAMake Headers
#include "math/xyz_vector.h"
#include "math/xyz_matrix.h"

int
are_floats_equal(
	float const,
	float const);

int
are_xyzVector_equal(
	Vector const &,
	Vector const &);

int
are_xyzVectors_equal(
	Vectors const &,
	Vectors const &
);

int
are_xyzMatrix_equal(
	Matrix const &,
	Matrix const &);

#endif /* defined(__REDESIGNC__Numeric_Test__) */
