//
//  structure.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__structure__
#define __RNAMake__structure__

#include <stdio.h>
#include "chain.h"
#include "types.h"
#include "transform.h"
#include "xyzMatrix.h"

class Structure {
public:
    Structure():
    chains_ ( Chains() ),
    dummy_ ( Point (0, 0, 0))
    {}
    
    Structure
    copy();
    
    ~Structure() {}
    //Structure

public:
    
    inline
    void
    _cache_coords() {
        coords_ = Points();
        org_coords_ = Points();
        for(auto const & a : atoms()) {
            coords_.push_back(a->coords());
            org_coords_.push_back(a->coords());
        }
        
    }
    
    void
    _build_chains(
        Residues &);
    
    inline
    void
    _update_coords(
        Points const & ncoords) {
        int i = 0;
        for (auto & a : atoms()) {
            a->coords(ncoords[i]);
            i++;
        }
    }
    
    inline
    void
    _restore_coords() {
        _update_coords(org_coords_);
    }
    
    inline
    void
    move(Point const & p) {
        for(auto & a : atoms()) {
            a->coords(a->coords() + p);
        }
    }
    
    inline
    void
    transform(Transform const & t) {
        Matrix r = t.rotation().transpose();
        Point trans = t.translation();
        for( auto & p : coords_) {
            dot_vector(r, p, dummy_);
            dummy_ += trans;
            p = dummy_;
        }
        
        
        _update_coords(coords_);
    }
    
    inline
    Beads
    get_beads(
        Residues const & excluded_res) {
        Beads beads;
        int found = 0;
        for ( auto const & r : residues()) {
            found = 0;
            for (auto const & er : excluded_res) {
                if( r == er) {found = 1; break;}
            }
            if (found) { continue; }
            for (auto const & b : r.get_beads()) {
                beads.push_back(b);
            }
        }
        return beads;
    }
    
    String
    to_pdb_str();
    
    String
    to_str();
    
    void
    to_pdb(String const);    
    
public: // getters
    inline
    Chains const &
    chains() { return chains_; }
    
    inline
    Residues const
    residues() {
        Residues residues;
        for (auto const & c : chains_) {
            for (auto const & r : c.residues()) {
                residues.push_back(r);
            }
        }
        return residues;
    }
    
    inline
    AtomOPs const
    atoms() {
        AtomOPs atoms;
        for (auto const & r : residues()) {
            for (auto const & a : r.atoms()) {
                if(a != NULL) {
                    atoms.push_back(a);
                }
            }
        }
        return atoms;
    }
    
    inline
    Residue &
    get_residue(
        int const & num,
        String const & chain_id,
        String const & i_code) {
        for( auto const & c : chains_) {
            for (auto & r : c.residues() ){
                if (num == r.num() && chain_id.compare(r.chain_id()) == 0 && i_code.compare(r.i_code()) == 0) {
                    return r;
                }
            }
        }
        throw "could not find residue";
    }
    
    inline
    Residue &
    get_residue(
        Uuid const & uuid) {
        for( auto const & c : chains_) {
            for (auto & r : c.residues() ){
                if ( r.uuid() == uuid) { return r; }
            }
        }
        throw "could not find residue";
    }

public: // setters
    
    inline
    void
    chains(Chains const & nchains) { chains_ = nchains; }
    
private:
    Chains chains_;
    Point dummy_; // resuable place in memory
    Points coords_;
    Points org_coords_;
    
};

Structure
str_to_structure(
    String const &,
    ResidueTypeSet const & rts);


#endif /* defined(__RNAMake__structure__) */
