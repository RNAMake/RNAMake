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
#include "uuid.h"
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
    copy() const {
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
        atoms_ ( AtomOPs() ),
        uuid_ ( Uuid() )
    {}
    
    Residue
    copy() const;
    
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
        float cutoff = 3.0) const {
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
    new_uuid() { uuid_ = Uuid(); }

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
    get_beads() const;

    String
    to_str() const;
    
    String
    to_pdb_str(int &) const;
    
    void
    to_pdb(String const);
 
    bool
    operator ==(const Residue& r) const {
        return uuid_ == r.uuid_;
    }

    
public: // setters
    inline
    void
    uuid(Uuid const & nuuid) { uuid_ = nuuid; }
    
    
public: // getters
    
    inline
    String const &
    name() const { return name_; }
    
    inline
    String const &
    chain_id() const { return chain_id_; }
    
    inline
    String const &
    i_code() const { return i_code_; }
    
    inline
    int const &
    num() const { return num_; }
    
    inline
    String const &
    short_name() const { return rtype_.short_name(); }
    
    inline
    AtomOPs const &
    atoms() const { return atoms_; }
    
    inline
    Uuid const &
    uuid() const { return uuid_; }
    
private:
    ResidueType rtype_;
    String name_, chain_id_, i_code_;
    int num_;
    AtomOPs atoms_;
    Uuid uuid_;
    
};


Residue
str_to_residue(
    String const & s,
    ResidueTypeSet const & rts);

typedef std::vector<Residue> Residues;
typedef std::shared_ptr<Residue> ResidueOP;
typedef std::vector<ResidueOP> ResidueOPs;




#endif /* defined(__RNAMake__residue__) */
