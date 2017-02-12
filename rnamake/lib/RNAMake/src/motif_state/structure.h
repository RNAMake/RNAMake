//
// Created by Joseph Yesselman on 1/29/17.
//

#ifndef RNAMAKE_MOTIF_STATE_STRUCTURE_H
#define RNAMAKE_MOTIF_STATE_STRUCTURE_H

#include "primitives/structure.h"
#include "motif_state/residue.h"
#include "motif_state/chain.h"

namespace state {

class Structure : public primitives::Structure<Chain, Residue> {
public:
    inline
    Structure(
            ResidueOPs const & res,
            Ints const & chain_cuts):
            primitives::Structure<Chain, Residue>(res, chain_cuts) {}

    inline
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
            String const & s):
            primitives::Structure<Chain, Residue>() {
        residues_ = ResidueOPs();
        Strings spl = split_str_by_delimiter(s, ";");
        for(int i = 0; i < spl.size()-1; i++) {
            residues_.push_back(std::make_shared<Residue>(spl[i]));
        }

        auto chain_cuts_spl = split_str_by_delimiter(spl.back(), " ");
        for(auto const & i : chain_cuts_spl) { chain_cuts_.push_back(std::stoi(i)); }
    }

    ~Structure() {}

public:

    String
    to_str() {
        String s;
        for (auto const & r : residues_) {
            s += r->to_str() + ";";
        }
        for (auto const & i : chain_cuts_) {
            s += std::to_string(i) + " ";
        }
        s += ";";
        return s;
    }

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
            Vector const & t,
            Point & dummy) {
        for(auto const & res : residues_) { res->fast_transform(r, t, dummy); }
    }

    ChainOPs
    get_chains() {
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

};

typedef std::shared_ptr<Structure> StructureOP;

}

#endif //RNAMAKE_MOTIF_STATE_STRUCTURE_H
