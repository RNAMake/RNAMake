//
//  motif_type.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "util/motif_type.h"

namespace util {

String const
type_to_str(MotifType const mtype) {
    if (mtype == MotifType::TWOWAY) { return "TWOWAY"; }
    else if (mtype == MotifType::NWAY) { return "NWAY"; }
    else if (mtype == MotifType::HAIRPIN) { return "HAIRPIN"; }
    else if (mtype == MotifType::T_T) { return "2X_TWOWAY"; }
    else if (mtype == MotifType::T_T_T) { return "3X_TWOWAY"; }
    else if (mtype == MotifType::TWOWAY_SEGMENTS) { return "TWOWAY_SEGMENTS"; }
    else if (mtype == MotifType::HELIX) { return "HELIX"; }
    else if (mtype == MotifType::UNKNOWN) { return "UNKNOWN"; }
    else if (mtype == MotifType::TCONTACT) { return "TCONTACT"; } // Added by CJ 08/20
    else { throw "cannot indentify type "; }
}

MotifType const
str_to_type(String const s) {
    if (s.compare("TWOWAY") == 0) { return MotifType::TWOWAY; }
    else if (s.compare("NWAY") == 0) { return MotifType::NWAY; }
    else if (s.compare("HAIRPIN") == 0) { return MotifType::HAIRPIN; }
    else if (s.compare("2X_TWOWAY") == 0) { return MotifType::T_T; }
    else if (s.compare("3X_TWOWAY") == 0) { return MotifType::T_T_T; }
    else if (s.compare("TWOWAY_SEGMENTS") == 0) { return MotifType::TWOWAY_SEGMENTS; }
    else if (s.compare("HELIX") == 0) { return MotifType::HELIX; }
    else if (s.compare("TCONTACT") == 0) { return MotifType::TCONTACT; } //added by CJ 08/20
    else if (s.compare("UNKNOWN") == 0) { return MotifType::UNKNOWN; }
    else { throw "cannot indentify str for type"; }

}

}



