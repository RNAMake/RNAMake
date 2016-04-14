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

/**
 * BeadType is an ENUM type. This is to specify which center of atoms each bead
 * represents.
 *
 * Phosphate (0):  P, OP1, OP2\n
 * Sugar (1):  O5',C5',C4',O4',C3',O3',C1',C2',O2'\n
 * Base  (2):  All remaining atoms
 */
enum BeadType { PHOS = 0, SUGAR = 1, BASE = 2 };


/**
 * Bead class stores information related to keeping track of steric clashes
 * between residues during building. They are never used outside the Residue class
 *
 */
class Bead {
public:
    /**
     * Empty constructor for Bead object.
     */
    Bead():
    center_(Point(0,0,0)),
    btype_(BeadType(0))
    {}
    
    /**
     * Standard constructor for Bead object.
     * @param   btype   type of bead (PHOS, SUGAR or BASE)
     * @param   center  the average 3D position of all atoms represented by bead
     */
    inline
    Bead(
        Point const & center,
        BeadType const btype):
        center_ ( center ),
        btype_ ( btype )
    {}
    
    /**
     * Copy constructor
     * @param   b   Bead object copying from
     */
    inline
    Bead(
        Bead const & b):
    center_(b.center_),
    btype_(b.btype_)
    {}
    
public: //accessors
    
    /**
     * Accessor for center_
     */
    inline
    Point
    center() const { return center_; }
    
    /**
     * Accessor for btype_
     */
    inline
    BeadType
    btype() const { return btype_; }

private:
    /**
     * private variable for the 3D coordinates of the center of atoms the bead represents
     */
    Point center_;
    
    /**
     * private variable of the type of the bead PHOS, SUGAR or BASE)
     */
    BeadType btype_;
    
};

typedef std::vector<Bead> Beads;

/**
 * Store residue information from pdb file, stores all Atom objects that
 * belong to residue. Implementation is designed to be extremely lightweight.
 *
 * To get an example instance of Residue please include:\n
 * #include "instances/structure_instances.hpp"
 *
 * Example of Common Usage:
 *
 * @code
 *  //grabbing example instance
 *  #include "instances/structure_instances.hpp"
 *  auto r = instances::residue();
 *  //simply getting name of residue, etc: A,G,C or U
 *  std::cout << r->name() << std::endl; //OUTPUT: "G"
 *
 *  //getting a specific atom from residue
 *  auto a = r->get_atom("C1'");
 *  std::cout << a->coords() << std::endl; //OUTPUT: "-23.806 -50.289  86.732"
 *
 *  auto beads = r->get_beads();
 *  std::cout << beads[0]->btype() << " " << beads[0]->center() << std::endl;
 *  //OUTPUT: "1 24.027 -48.5001111111 86.368"
 * 
 *  //get PDB formatted String (can also write to file with r->to_pdb("test.pdb") )
 *  std::cout << r->to_pdb_str() << std::endl; //OUTPUT -->
 *  ATOM      1 O5'  G   A 103     -26.469 -47.756  84.669  1.00  0.00
 *  ATOM      2 C5'  G   A 103     -25.050 -47.579  84.775  1.00  0.00
 *  ATOM      3 C4'  G   A 103     -24.521 -48.156  86.068  1.00  0.00
 *  ATOM      4 O4'  G   A 103     -24.861 -49.568  86.118  1.00  0.00
 *  ATOM      5 C3'  G   A 103     -23.009 -48.119  86.281  1.00  0.00
 *  ATOM      6 O3'  G   A 103     -22.548 -46.872  86.808  1.00  0.00
 *  .
 *  .
 *  .
 * @endcode
 */

class Residue {
public:
    
    /**
     * Standard constructor for Bead object.
     * @param   rtype   residue type, dictates residue topology and information
     * @param   name    name of residue i.e. "GUA", "CYT", etc
     * @param   num     residue num 
     * @param   chain_id    the id of the chain i.e. "A", "B", only one character
     * @param   i_code  residue insertion code, usually nothing ("")
     */
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
    
    /**
     * Copy constructor
     * @param   r   Residue object copying from
     */
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
    
    /**
     * Construction from String, used in reading data from files
     * @param   s   string generated from to_str()
     * @see to_str()
     */
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

    
    /**
     * Empty deconstructor
     */
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
    
    Beads
    get_beads() const;

    String
    to_str() const;
    
    inline
    String
    to_pdb_str(int & acount) {
        return to_pdb_str(acount, -1, "");
    }
    
    String
    to_pdb_str(
        int &,
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
    
    inline
    void
    uuid(Uuid const & uuid) { uuid_ = uuid; }
    
    
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
