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
#include "motif_state/structure.h"
#include "structure/chain.fwd.h"
#include "structure/chain.h"
#include "structure/residue.h"
#include "structure/pdb_parser.h"

class Structure : public primitives::Structure<Chain, Residue> {
public:

    Structure(
            ResidueOPs const & res,
            Ints const & chain_cuts):
            primitives::Structure<Chain, Residue>(res, chain_cuts) {}

    Structure(
            Structure const & s,
            int new_uuid = 0):
            primitives::Structure<Chain, Residue>() {

        residues_ = ResidueOPs();
        for(auto const & r : s.residues_) {
            residues_.push_back(std::make_shared<Residue>(*r, new_uuid));
        }
        chain_cuts_ = s.chain_cuts_;

    }

    Structure(
            String const & s,
            ResidueTypeSet const & rts):
            primitives::Structure<Chain, Residue>() {
        residues_ = ResidueOPs();
        Strings spl = split_str_by_delimiter(s, ";");
        for(int i = 0; i < spl.size()-1; i++) {
            residues_.push_back(std::make_shared<Residue>(spl[i], rts));
        }

        auto chain_cuts_spl = split_str_by_delimiter(spl.back(), " ");
        for(auto const & i : chain_cuts_spl) { chain_cuts_.push_back(std::stoi(i)); }

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

    ChainOPs
    get_chains() const {
        auto pos = 0;
        auto res = ResidueOPs();
        auto chains = ChainOPs();
        auto i = -1;
        for(auto const & r : residues_) {
            i++;
            if(chain_cuts_[pos] == i) {
                auto c = std::make_shared<Chain>(res);
                chains.push_back(c);
                res = ResidueOPs();
                res.push_back(r);
                pos++;
            }
            else {
                res.push_back(r);
            }
        }

        if(res.size() > 0) {
            chains.push_back(std::make_shared<Chain>(res));
        }

        return chains;
    }

    state::StructureOP
    get_state() {
        auto residues = state::ResidueOPs();
        for(auto const & r : residues_) {
            residues.push_back(r->get_state());
        }
        return std::make_shared<state::Structure>(residues, chain_cuts_);
    }


};

typedef std::shared_ptr<Structure> StructureOP;

bool
are_structures_equal(
        StructureOP const & s1,
        StructureOP const & s2,
        int check_uuids = 1);

StructureOP
structure_from_pdb(
        String const &,
        ResidueTypeSet const &);


#endif /* defined(__RNAMake__structure__) */
