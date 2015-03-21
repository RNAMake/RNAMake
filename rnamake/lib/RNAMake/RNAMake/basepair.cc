//
//  basepair.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "basepair.h"


Basepair
Basepair::copy() {
    Basepair cbp ( res1_, res2_, bp_state_.r(), bp_type_);
    cbp.flipped_ = flipped_;
    cbp.uuid_ = uuid_;
    return cbp;
}

String const
Basepair::to_str() const {
    std::stringstream ss;
    ss << name() << "," << bp_state_.to_str() << "," << bp_type_ << ",0," << flipped_;
    return ss.str();
}

String const
Basepair::to_pdb_str() const {
    String s;
    int acount = 1;
    for (auto const & r : residues()) {
        r->to_pdb_str(acount);
    }
    return s;
}

void
Basepair::to_pdb(String const fname) const {
    std::ofstream out;
    out.open(fname.c_str());
    String s = to_pdb_str();
    out << s << std::endl;
    out.close();
}

bool
wc_bp(BasepairOP const & bp) {
    String bp_str = bp->res1()->short_name() + bp->res2()->short_name();
    if(bp_str.compare("GC") == 0) { return true; }
    if(bp_str.compare("CG") == 0) { return true; }
    if(bp_str.compare("AU") == 0) { return true; }
    if(bp_str.compare("UA") == 0) { return true; }
    return false;
}

bool
gu_bp(BasepairOP const & bp) {
    String bp_str = bp->res1()->short_name() + bp->res2()->short_name();
    if(bp_str.compare("GU") == 0) { return true; }
    if(bp_str.compare("UG") == 0) { return true; }
    return false;
}

