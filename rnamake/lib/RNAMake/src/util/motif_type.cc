//
//  motif_type.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "util/motif_type.h"

String const
type_to_str(MotifType const mtype) {
    if     (mtype == TWOWAY             ) {   return "TWOWAY";          }
    else if(mtype == NWAY               ) {   return "NWAY";            }
    else if(mtype == HAIRPIN            ) {   return "HAIRPIN";         }
    else if(mtype == T_T                ) {   return "2X_TWOWAY";       }
    else if(mtype == T_T_T              ) {   return "3X_TWOWAY";       }
    else if(mtype == TWOWAY_SEGMENTS    ) {   return "TWOWAY_SEGMENTS"; }
    else if(mtype == HELIX              ) {   return "HELIX";           }
    else if(mtype == UNKNOWN            ) {   return "UNKNOWN";         }
    else { throw "cannot indentify type "; }
}

MotifType const
str_to_type(String const s) {
    if     (s.compare("TWOWAY") == 0           ) {   return TWOWAY;          }
    else if(s.compare("NWAY") == 0             ) {   return NWAY;            }
    else if(s.compare("HAIRPIN") == 0          ) {   return HAIRPIN;         }
    else if(s.compare("2X_TWOWAY") == 0        ) {   return T_T;             }
    else if(s.compare("3X_TWOWAY") == 0        ) {   return T_T_T;           }
    else if(s.compare("TWOWAY_SEGMENTS") == 0  ) {   return TWOWAY_SEGMENTS; }
    else if(s.compare("HELIX") == 0            ) {   return HELIX;           }
    else if(s.compare("UNKNOWN") == 0          ) {   return UNKNOWN;         }
    else { throw "cannot indentify str for type"; }

}




