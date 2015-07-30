/**
 * @file
 * @author  Joseph D. Yesselman
 * @version 1.0
 *
 * @section DESCRIPTION
 *
 *  stores atomic information from pdb file, design is to be extremely
 *  lightweight only storing the atom name and coordinates.
 */

#ifndef __RNAMake__atom__
#define __RNAMake__atom__

#include <stdio.h>

//RNAMake Headers
#include "base/types.h"
#include "math/xyz_vector.h"

class Atom {
public:
    friend class Structure;
    
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
    Point const
    coords() const { return coords_; }

    inline
    Point const &
    coords() { return coords_; }

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

typedef std::shared_ptr<Atom> AtomOP;
typedef std::vector<AtomOP> AtomOPs;

typedef std::shared_ptr<const Atom> AtomCOP;
typedef std::vector<AtomCOP>  AtomCOPs;




#endif /* defined(__RNAMake__atom__) */
