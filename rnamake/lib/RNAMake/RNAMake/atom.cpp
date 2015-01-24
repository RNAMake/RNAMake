//
//  atom.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/23/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#include "atom.h"
#include "FileIO.h"
#include "xyzVector.h"

Atom::Atom(
    String const & name,
    Point const & coords):
    name_ ( name ),
    coords_ ( coords )
{}

Atom
Atom::copy() {
    return Atom(name_, coords_);
}

String
Atom::to_str() {
    return name_ + " " + vector_to_str(coords_);
}

String
Atom::to_pdb_str(
    int acount) {
    
    char buffer [200];
    std::sprintf(buffer, "ATOM %6d  P   C   A   1 %11.3f%8.3f%8.3f  1.00 62.18           P\n", acount, coords_[0], coords_[1], coords_[2]);    
    return String(buffer);
    
}

Atom
str_to_atom(
    String const & s) {
    
    Strings spl = split_str_by_delimiter(s, " ");
    Point coords(std::stof(spl[1]), std::stof(spl[2]), std::stof(spl[3]));
    return Atom(spl[0], coords);
}


