//
//  residue.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__residue__
#define __RNAMake__residue__

#include <stdio.h>
#include <uuid/uuid.h>
#include "residue_type.h"
#include "residue_type_set.h"
#include "types.h"
#include "atom.h"

enum BeadType { PHOS = 0, SUGAR = 1, BASE = 2 };

class Bead {
public:
    Bead() {}
    
    inline
    Bead(
        Point const & center,
        BeadType const btype):
        center_ ( center ),
        btype_ ( btype )
    {}
    
    inline
    Bead
    copy() {
        return Bead(center_, btype_);
    }

public:
    
    inline
    Point
    center() { return center_; }
    
    inline
    BeadType
    btype() { return btype_; }

private:
    Point center_;
    BeadType btype_;
    
};

class Residue {
public:
    Residue() {}
    
    inline
    Residue(
        ResidueType const & rtype,
        String const & name,
        int const & num,
        String const & chain_id,
        String const & i_code):
        rtype_ ( rtype ),
        name_ ( name ),
        num_ ( num ),
        chain_id_ ( chain_id ),
        i_code_ ( i_code ),
        atoms_ ( AtomOPs() )
    {  uuid_generate_random(uuid_); }
    
    ~Residue() {}

public:
    
    void
    setup_atoms(
        AtomOPs &);
    
public:
    
    inline
    AtomOP
    const &
    get_atom(
        String const & name) {
        int index = rtype_.atom_pos_by_name(name);
        if( index == -1) {
            throw("cannot find atom name " + name);
        }
        return atoms_[index];
    }
    
    
private:
    ResidueType rtype_;
    String name_, chain_id_, i_code_;
    int num_;
    uuid_t uuid_;
    AtomOPs atoms_;
    
};

inline
Residue
str_to_residue(
    String const & s,
    ResidueTypeSet const & rts) {
    Strings spl = split_str_by_delimiter(s, ",");
    ResidueType rtype = rts.get_rtype_by_resname(spl[0]);
    Residue r(rtype, spl[1], std::stoi(spl[2].c_str()), spl[3], spl[4]);
    AtomOPs atoms;
    int i = 5;
    while ( i < spl.size()-1 ) {
        if( spl[i].length() == 1) {
            atoms.push_back( NULL );
        }
        else {
            AtomOP aop ( new Atom ( str_to_atom(spl[i]) ));
            atoms.push_back(aop);
        }
        i++;
    }
    r.setup_atoms(atoms);
    return r;
}

typedef std::vector<Residue> Residues;





#endif /* defined(__RNAMake__residue__) */
