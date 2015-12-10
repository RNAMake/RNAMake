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

//RNAMake Headers
#include "base/types.h"
#include "math/transform.h"
#include "math/xyz_matrix.h"
#include "structure/chain.fwd.h"
#include "structure/chain.h"
#include "structure/residue.h"
#include "structure/pdb_parser.h"

class Structure {
public:
    inline
    Structure():
    chains_ (ChainOPs())
    {}
    
    inline
    Structure(
        ChainOPs const & chains):
    chains_ (chains)
    {}
    
    Structure(
        String const & path) {
        PDBParser pdb_parser;
        auto residues = pdb_parser.parse(path);
        chains_ = ChainOPs();
        connect_residues_into_chains(residues, chains_);
    }
    
    Structure(
        Structure const & s) {
        chains_ = ChainOPs(s.chains_.size());
        int i = 0;
        for (auto const & c : s.chains_) {
            chains_[i] = std::make_shared<Chain>(*c);
            i++;
        }
    }
    
    Structure(
        String const & s,
        ResidueTypeSet const & rts) {
        chains_ = ChainOPs();
        Strings spl = split_str_by_delimiter(s, ":");
        for( auto const & c_str : spl) {
            chains_.push_back(std::make_shared<Chain>(c_str, rts));
        }
    }
    
    ~Structure() {}
    
public:
    
    void
    renumber();
    
    inline
    Beads
    get_beads(
        ResidueOPs const & excluded_res) {
        Beads beads;
        int found = 0;
        for ( auto const & r : residues()) {
            found = 0;
            for (auto const & er : excluded_res) {
                if( r->uuid() == er->uuid()) {found = 1; break;}
            }
            if (found) { continue; }
            for (auto const & b : r->get_beads()) {
                beads.push_back(b);
            }
        }
        return beads;
    }
    
    ResidueOP const
    get_residue(
        int const & ,
        String const & ,
        String const & );
    
    ResidueOP const
    get_residue(
        Uuid const &);

    ResidueOPs const
    residues();
    
    inline
    AtomOPs const
    atoms() {
        AtomOPs atoms;
        for (auto const & r : residues()) {
            for (auto const & a : r->atoms()) {
                if(a.get() != NULL) {
                    atoms.push_back(a);
                }
            }
        }
        return atoms;
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
        for( auto & a : atoms() ) {
            dot_vector(r, a->coords(), dummy_);
            dummy_ += trans;
            a->coords(dummy_);
        }
    }
    
      String
    to_pdb_str();
    
    String
    to_str();
    
    void
    to_pdb(String const);
    
public: // getters
    
    inline
    ChainOPs const &
    chains() { return chains_; }

    
private:
    ChainOPs chains_;
    Point dummy_; // resuable place in memory
    Points coords_;
    Points org_coords_;
    
};

typedef std::shared_ptr<Structure> StructureOP;


#endif /* defined(__RNAMake__structure__) */
