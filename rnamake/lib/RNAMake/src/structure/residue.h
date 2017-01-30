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
#include "util/bead.h"
#include "primitives/residue.h"
#include "motif_state/residue.h"
#include "structure/atom.h"
#include "structure/residue_type.h"
#include "structure/residue_type_set.h"

Point
calc_center(AtomOPs const &);

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

class Residue : public primitives::Residue {
public:

    /**
     * Standard constructor for Residue object.
     * @param   rtype   residue type, dictates residue topology and information
     * @param   name    name of residue i.e. "GUA", "CYT", etc
     * @param   num     residue num 
     * @param   chain_id    the id of the chain i.e. "A", "B", only one character
     * @param   i_code  residue insertion code, usually nothing ("")
     */
    Residue(
            AtomOPs const &,
            ResidueTypeOP const &,
            String const &,
            int const &,
            String const &,
            String const &);

    /**
     * Copy constructor
     * @param   r   Residue object copying from
     */
    Residue(
            Residue const & r,
            int new_uuid = 0,
            int build_beads = 1);

    Residue(
            Residue const & r,
            Uuid const & given_uuid,
            int build_beads = 1);

    /**
     * Construction from String, used in reading data from files
     * @param   s   string generated from to_str()
     * @see to_str()
     */
    Residue(
            String const &,
            ResidueTypeSet const &);

    /**
     * Empty deconstructor
     */
    ~Residue() {}

public: //iterator
    typedef typename AtomOPs::iterator iterator;
    typedef typename AtomOPs::const_iterator const_iterator;

    iterator begin() { return atoms_.begin(); }
    iterator end()   { return atoms_.end(); }

    const_iterator begin() const { return atoms_.begin(); }
    const_iterator end()   const { return atoms_.end(); }

public:
    /**
     * put atoms in correct positon in internal atom list, also corrects some
     * named atom names to their correct name
     * @param   atoms    list of atom objects that are to be part of this residue
     */
    void
    setup_atoms(
            AtomOPs const &);

public:

    /**
     * get atom object by its name
     * @param   name   name of atom
     *
     * examples:
     * @code 
     *  auto a = r->get_atom("C1'");
     *  std::cout << a->coords() << std::endl;
     *  //OUTPUT -23.806 -50.289  86.732
     * @endcode
     */
    inline
    AtomOP const &
    get_atom(String const & name) const {
        auto result = rtype_->is_valid_atom(name);
        if(!result) {
            throw ResidueException("cannot find atom name " + name);
        }
        auto index = rtype_->atom_index(name);
        return get_atom(index);
    }

    inline
    AtomOP const &
    get_atom(int index) const{
        if(index >= atoms_.size()) {
            throw ResidueException("cannot get atom with index: " + std::to_string(index));
        }
        if(atoms_[index] == nullptr) {
            throw ResidueException("cannot get atom wiht index: " + std::to_string(index) +
                                   " it is not initialized");
        }

        return atoms_[index];
    }

    inline
    bool
    has_atom(String const & name) const {
        auto result = rtype_->is_valid_atom(name);
        if(!result) { return false; }

        auto index = rtype_->atom_index(name);
        return has_atom(index);
    }

    inline
    bool
    has_atom(int index) const {
        if(index >= atoms_.size()) { return false; }
        if(atoms_[index] == nullptr) { return false; }

        return true;
    }

    /**
     * Determine if another residue is connected to this residue, returns 0
     * if res is not connected to self, returns 1 if connection is going
     *  from 5' to 3' and returns -1 if connection is going from 3' to 5'
     * @param   res another residue
     * @param   cutoff  distance to be considered connected, default: 3 Angstroms
     */
    inline
    int
    connected_to(
            Residue const & res,
            float cutoff = 3.0) const {

        // 5' to 3'
        if(has_atom("O3'") && res.has_atom("P")) {
            auto o3_atom = get_atom("O3'");
            auto p_atom  = res.get_atom("P");
            if (o3_atom != nullptr && p_atom != nullptr) {
                if (o3_atom->coords().distance(p_atom->coords()) < cutoff) {
                    return 1;
                }
            }
        }

        // 3' to 5'
        if(has_atom("P") && res.has_atom("O3'")) {
            auto o3_atom = res.get_atom("O3'");
            auto p_atom = get_atom("P");
            if (o3_atom != nullptr && p_atom != nullptr) {
                if (o3_atom->coords().distance(p_atom->coords()) < cutoff) {
                    return -1;
                }
            }
        }

        return 0;
    }

    /**
     * wrapper for to_pdb_str(int &, int, String const &) when one does not care about
     * renumbering atoms and residue
     */
    inline
    String
    to_pdb_str(int & acount) {
        return to_pdb_str(acount, -1, "");
    }

    inline
    void
    build_beads() { beads_ = _get_beads(); }

    inline
    void
    remove_beads() { beads_ = Beads(); }

    inline
    void
    move(Point const & p) {
        for(auto & a : atoms_) {
            if(a != nullptr) { a->move(p); }
        }
        for(auto & b : beads_) { b.move(p); }
    }

    inline
    void
    transform(Transform const & t) {
        for(auto & a : atoms_) {
            if(a != nullptr) { a->transform(t); }
        }
        for(auto & b : beads_) { b.transform(t); }
    }

    inline
    void
    fast_transform(
            Matrix const & r,
            Vector const & t) {
        for(auto & a : atoms_) {
            if(a != nullptr) { a->fast_transform(r, t, rtype_->dummy_coords()); }
        }

        //for(auto & b : beads_) { b.fast_transform(r, t); }
    }

public:

    /**
     * stringifes residue object
     * @code
     *  #include "instances/structure_instances.hpp"
     *  auto r = instances::residue();
     *  std::cout << r->to_str() << std::endl;
     *  //OUTPUT
     *  GUA,G,103,A,,N,N,N,O5' -26.469 -47.756 84.669,C5' -25.05 -47.579 84.775,C4' -24.521 -48.156
     *  86.068,O4' -24.861 -49.568 86.118,C3' -23.009 -48.119 86.281,O3' -22.548 -46.872 86.808,C1'
     *  -23.806 -50.289 86.732,C2' -22.812 -49.259 87.269,O2' -23.167 -48.903 88.592,N1 -19.538 -52.485
     *  85.025,C2 -19.717 -51.643 86.097,N2 -18.624 -51.354 86.809,N3 -20.884 -51.124 86.445,C4 -21.881
     *  -51.521 85.623,C5 -21.811 -52.356 84.527,C6 -20.546 -52.91 84.164,O6 -20.273 -53.677 83.228,N7
     *  -23.063 -52.513 83.947,C8 -23.858 -51.786 84.686,N9 -23.21 -51.159 85.722,
     * @endcode
     */
    String
    to_str() const;

    /**
   * writes a PDB string formmated verision of this Residue object to file
   * @param  filename of output PDB file
   */
    void
    to_pdb(String const);

    /**
 * returns pdb formatted string of residue's coordinate information
 * @param   acount  current atom index, default: 1
 * @param   rnum    starting residue number
 * @param   chain_id    the chain id of the chain, i.e. "A", "B" etc
 * @code
 *  #include "instances/structure_instances.hpp"
 *  auto r = instances::residue();
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
    String
    to_pdb_str(
            int &,
            int,
            String const &) const;

    state::ResidueOP
    get_state();

public: // getters

    /**
     * getter for the one letter residue type
     */
    inline
    String
    short_name() const { return rtype_->short_name(); }

    inline
    Beads const &
    beads() { return beads_; }

    inline
    size_t
    num_beads() { return beads_.size(); }

    inline
    Point
    center() { return calc_center(atoms_); }

private:
    /**
 * Generates steric beads required for checking for steric clashes between
 * motifs. Each residues has three beads modeled after the typical three
 * bead models used in coarse grain modeling. The three beads are:
 *
 * @code
 *  #include "instances/structure_instances.hpp"
 *  auto r = instances::residue();
 *  auto beads = r->get_beads();
 *  //Test instance is the first residue in a chain with no phosphate so it has
 *  //only two beads
 *  std::cout << beads.size() << std::endl;
 *  //OUTPUT 2
 *  std::cout << beads[0].btype() == BeadType::SUGAR << std::endl;
 *  //OUTPUT 1
 * @endcode
 */
    Beads
    _get_beads() const;

private:
    /**
     * residue type object which explains what atoms in belong in this residue.
     */
    ResidueTypeOP rtype_;

    /**
     * vector of the atom objects that belong to this residue
     */
    AtomOPs atoms_;

    /**
    * vector of the beads represents sterics for Phosphate, Sugars and Base
    */
    Beads beads_;

};

/**
 * Shared pointer typedef for Residue. Only use shared pointers!
 */
typedef std::shared_ptr<Residue> ResidueOP;

/**
 * Typedef of a vector of shared pointer vectors, only used this.
 */
typedef std::vector<ResidueOP> ResidueOPs;


#endif /* defined(__RNAMake__residue__) */
