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

typedef std::vector<Bead> Beads;

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
    
    Residue
    copy();
    
    ~Residue() {}

public:
    
    void
    setup_atoms(
        AtomOPs &);
    
public:
    
    inline
    AtomOP
    const &
    get_atom (
        String const & name) const {
        int index = rtype_.atom_pos_by_name(name);
        if( index == -1) {
            throw("cannot find atom name " + name);
        }
        return atoms_[index];
    }
    
    inline
    int
    connected_to(
        Residue const & res,
        float cutoff) const {
        String o3 = "O3'", p = "P";
        
        // 5' to 3'
        AtomOP o3_atom = get_atom(o3), p_atom = res.get_atom(p);
        if(o3_atom != NULL && p_atom != NULL) {
            if( o3_atom->coords().distance(p_atom->coords()) < cutoff) {
                return 1;
            }
        }
        
        // 3' to 5'
        o3_atom = res.get_atom(o3); p_atom = get_atom(p);
        if(o3_atom != NULL && p_atom != NULL) {
            if( o3_atom->coords().distance(p_atom->coords()) < cutoff) {
                return -1;
            }
        }
        
        return 0;
    }
    
    void
    new_uuid() { uuid_generate_random(uuid_); }

    inline
    int
    acount() {
        int count = 0;
        for (auto const & a : atoms_) {
            if (a != NULL) { count++; }
        }
        return count;
    }
    
    Beads
    get_beads();

    String
    to_str();
    
    String
    to_pdb_str(int &);
    
    void
    to_pdb(String const);
    
public: // setters
    inline
    void
    uuid(uuid_t const & nuuid) { uuid_copy(uuid_, nuuid); }
    
public: // getters
    
    inline
    String const &
    name() { return name_; }
    
    inline
    String const &
    chain_id() { return chain_id_; }
    
    inline
    String const &
    i_code() { return i_code_; }
    
    inline
    int const &
    num() { return num_; }
    
    inline
    String const &
    short_name() const { return rtype_.short_name(); }
    

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
    while ( i < spl.size() ) {
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
