//
//  atom.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/23/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

// RNAMake Headers
#include <math/vector_3.hpp>
#include <structure/all_atom/atom.h>

namespace structure::all_atom {

String Atom::get_str() const {
  return _name + " " + _coords.get_str();
}

} // namespace structure::all_atom
