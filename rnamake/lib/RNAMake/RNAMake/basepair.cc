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