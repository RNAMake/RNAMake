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
    
    inline
    Atom(
        String const & name,
        Point const & coords):
        name_ ( name ),
        coords_ ( coords )
    {}
    
    inline
    Atom
    copy() const {
        return Atom(name_, coords_);
    }
    
    String to_str();
    String to_pdb_str(int);
public: //accessors
    
    inline
    String const &
    name() { return name_; }
    
    inline
    String const &
    name() const { return name_; }
    
    inline
    Point const &
    coords() { return coords_; }
    
    inline
    Point const
    coords() const { return coords_; }
public: // setters
    
    inline
    void
    coords(
        Point const & ncoords) {
        coords_ = ncoords;
    }
    
    inline
    void
    name(
        String const & nname) {
        name_ = nname;
    }
    
    
private:
    String name_;
    Point coords_;
    
};

Atom
str_to_atom(
    String const &);

typedef std::vector<Atom> Atoms;
typedef std::shared_ptr<Atom> AtomOP;
typedef std::vector<AtomOP> AtomOPs;

#endif /* defined(__RNAMake__atom__) */
