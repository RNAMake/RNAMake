//
//  atom.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/23/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#include "structure/atom.h"

//RNAMake Headers
#include "math/xyz_vector.h"

namespace structure {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// output functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

String
Atom::to_str() {
    return _name + " " + vector_to_str(_coords);
}

String
Atom::to_pdb_str(
        int acount) {

    char buffer[200];
    std::sprintf(buffer, "ATOM %6d  P   C   A   1 %11.3f%8.3f%8.3f  1.00 62.18           P\n", acount, _coords[0],
                 _coords[1], _coords[2]);
    return String(buffer);

}

}
