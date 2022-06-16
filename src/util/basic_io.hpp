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

#include "math/vector_3.hpp"

namespace util {

String points_to_pdb_str(math::Vector3s const &);

void points_to_pdb(String const &, math::Vector3s const &);

} // namespace util

#endif /* basic_io_hpp */
