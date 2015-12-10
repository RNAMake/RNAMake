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

//RNAMake Headers
#include "base/types.h"
#include "util/uuid.h"
#include "structure/atom.h"
#include "structure/residue_type.h"
#include "structure/residue_type_set.h"

Point
center(AtomOPs const &);

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
    Bead(
        Bead const & b):
    center_(b.center_),
    btype_(b.btype_)
    {}
    
public:
    
    inline
    Point
    center() const { return center_; }
    
    inline
    BeadType
    btype() const { return btype_; }

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
    
    Residue(
        Residue const & r):
    rtype_(r.rtype_),
    name_(r.name_),
    num_(r.num_),
    chain_id_(r.chain_id_),
    i_code_(r.i_code_),
    atoms_ ( AtomOPs(r.atoms().size()) ),
    uuid_(r.uuid_) {
        int i = -1;
        for(auto const & a : r.atoms_) {
            i++;
            if(a  == nullptr) { continue; }
            atoms_[i] = std::make_shared<Atom>(*a);
        }
    }
    
    Residue(
        String const & s,
        ResidueTypeSet const & rts) {
        
        Strings spl = split_str_by_delimiter(s, ",");
        rtype_      = rts.get_rtype_by_resname(spl[0]);
        name_       = spl[1];
        num_        = std::stoi(spl[2]);
        chain_id_   = spl[3];
        i_code_     = spl[4];
        atoms_      = AtomOPs();
        auto atoms  = AtomOPs();
        int i = 5;
        while ( i < spl.size() ) {
            if( spl[i].length() == 1) {
                atoms.push_back( nullptr );
            }
            else {
                atoms.push_back(std::make_shared<Atom>(spl[i]));
            }
            i++;
        }        
        setup_atoms(atoms);
    }

    
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
        if(o3_atom != nullptr && p_atom != nullptr) {
            if( o3_atom->coords().distance(p_atom->coords()) < cutoff) {
                return 1;
            }
        }
        
        // 3' to 5'
        o3_atom = res.get_atom(o3); p_atom = get_atom(p);
        if(o3_atom != nullptr && p_atom != nullptr) {
            if( o3_atom->coords().distance(p_atom->coords()) < cutoff) {
                return -1;
            }
        }
        
        return 0;
    }
    
    void
    new_uuid() {
        uuid_ = Uuid();
    }

    inline
    int
    acount() {
        int count = 0;
        for (auto const & a : atoms_) {
            if (a.get() != NULL) { count++; }
        }
        return count;
    }
    
    Beads
    get_beads() const;

    String
    to_str() const;
    
    String
    to_pdb_str(int & acount) {
        return to_pdb_str(acount, -1, "");
    }
    
    String
    to_pdb_str(int &,
               int,
               String const &) const;
    
    void
    to_pdb(String const);
 
    bool
    operator ==(const Residue& r) const {
        return uuid_ == r.uuid_;
    }

    
public: // setters
    inline
    void
    num(int nnum) { num_ = nnum; }
    
    inline
    void
    chain_id(String const & nchain_id) { chain_id_ = nchain_id; }
    
    
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


typedef std::shared_ptr<Residue> ResidueOP;
typedef std::vector<ResidueOP>   ResidueOPs;




#endif /* defined(__RNAMake__residue__) */
