//
//  basic_io.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 3/1/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef basic_io_hpp
#define basic_io_hpp

#include <stdio.h>

#include "math/xyz_vector.h"

String
points_to_pdb_str(
    math::Points const &);

void
points_to_pdb(
    String const &,
    math::Points const &);

#endif /* basic_io_hpp */
