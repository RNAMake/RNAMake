//
//  atom.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/23/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__atom__
#define __RNAMake__atom__

#include <stdio.h>
#include "types.h"
#include "xyzVector.h"

class Atom {
public:
    Atom() {}
    Atom(String const &, Point const &);
public:
    Atom copy();
    String to_str();
    String to_pdb_str(int);
public: //accessors
    String const &
    name() { return name_; }
    
    String const &
    name() const { return name_; }
    
    Point const &
    coords() { return coords_; }
    
    Point const
    coords() const { return coords_; }
public: // setters
    void
    coords(
        Point const & ncoords) {
        coords_ = ncoords;
    }
    
private:
    String name_;
    Point coords_;
    
};

#endif /* defined(__RNAMake__atom__) */
