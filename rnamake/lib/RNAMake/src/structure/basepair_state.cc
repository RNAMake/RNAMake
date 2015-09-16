//
//  BasepairState.cpp
//  REDESIGNC
//
//  Created by Joseph Yesselman on 9/28/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#include "util/settings.h"
#include "util/file_io.h"
#include "structure/basepair_state.h"

BasepairState
str_to_basepairstate(
	String const & s) {
	
	//always d,r,sugars order most be kept!!
	
	Strings strs = split_str_by_delimiter(s, ";");
	if(strs.size() < 3) {
		throw "cannot load BasepairState from String, not the right number of elements\n";
	}
	
	Vector d = vector_from_str(strs[0]);
	Matrix r = matrix_from_str(strs[1]);
	Vectors sug = vectors_from_str(strs[2]);
	
	BasepairState bp(d,r,sug);
	return bp;
}

std::ostream&
operator <<(
	std::ostream& stream,
	const BasepairState& bpstate) {
	
	stream << "<BasepairState Object(\n";
	stream << "\tOrigin\n";
	stream << bpstate.d() << "\n";
	stream << "\tRotation\n";
	stream << bpstate.r() << "\n";
	stream << "\tSugars\n";
	stream << bpstate.sugars()[0] << "\n";
	stream << bpstate.sugars()[1] << "\n";
	return stream;
	
}

BasepairState
get_ref_bp_state() {
	String path = resources_path() + "/ref_bp_state.dat";
	std::ifstream input;
	String line;
	
	input.open(path);
	getline(input,line);
	input.close();
	
	return str_to_basepairstate(line);
	
}

float
get_bpstate_rotation_diff(
	BasepairState const & bp1,
	BasepairState const & bp2) {
	
	float r_diff   = bp1.r().difference(bp2.r());
	Matrix flipped = bp2.r().get_flip_orientation();
	float r_diff_2 = bp1.r().difference(flipped);

	if(r_diff > r_diff_2) {
		r_diff = r_diff_2;
	}
	
	return r_diff;
	
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

