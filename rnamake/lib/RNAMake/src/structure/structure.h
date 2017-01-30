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
#include "primitives/structure.h"
#include "structure/chain.fwd.h"
#include "structure/chain.h"
#include "structure/residue.h"
#include "structure/pdb_parser.h"

class Structure : public primitives::Structure<Chain, Residue> {
public:

    Structure(
            ChainOPs const & chains):
            primitives::Structure<Chain, Residue>(chains) {}

    Structure(
            Structure const & s,
            int new_uuid = 0):
            primitives::Structure<Chain, Residue>() {

        chains_ = ChainOPs(s.chains_.size());
        int i = 0;
        for (auto const & c : s.chains_) {
            chains_[i] = std::make_shared<Chain>(*c, new_uuid);
            i++;
        }

        for(auto const & c : chains_) {
            for(auto const & r : *c) { residues_.push_back(r); }
        }
    }

    Structure(
            String const & s,
            ResidueTypeSet const & rts):
            primitives::Structure<Chain, Residue>() {
        chains_ = ChainOPs();
        Strings spl = split_str_by_delimiter(s, ":");
        for (auto const & c_str : spl) {
            chains_.push_back(std::make_shared<Chain>(c_str, rts));
        }

        for(auto const & c : chains_) {
            for(auto const & r : *c) { residues_.push_back(r); }
        }
    }

    ~Structure() {}

public:

    inline
    void
    move(Point const & p) {
        for(auto const & r : residues_) { r->move(p); }
    }

    inline
    void
    transform(Transform const & t) {
       for(auto const & r : residues_) { r->transform(t); }
    }

    inline
    void
    fast_transform(
            Matrix const & r,
            Point const & t) {
        for(auto const & res : residues_) { res->fast_transform(r, t); }
    }

    String
    to_pdb_str(
            int renumber = -1);

    String
    to_str();

    void
    to_pdb(
            String const,
            int renumber = -1);

public: // getters



};

typedef std::shared_ptr<Structure> StructureOP;

StructureOP
structure_from_pdb(
        String const &,
        ResidueTypeSet const &);


#endif /* defined(__RNAMake__structure__) */
