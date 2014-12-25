//
//  atom.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/23/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#include "atom.h"

Atom::Atom(
    String const & name,
    Point const & coords):
    name_ ( name ),
    coords_ ( coords )
{}