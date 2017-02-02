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
    Structure(ChainOPs const & chains):
            primitives::Structure<Chain, Residue>(chains) {}

    inline
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

    inline
    Structure(
            String const & s):
            primitives::Structure<Chain, Residue>() {
        chains_ = ChainOPs();
        Strings spl = split_str_by_delimiter(s, ":");
        for (auto const & c_str : spl) {
            chains_.push_back(std::make_shared<Chain>(c_str));
        }

        for(auto const & c : chains_) {
            for(auto const & r : *c) { residues_.push_back(r); }
        }
    }

    ~Structure() {}

public:

    String
    to_str() {
        String s;
        for (auto const & c : chains_) {
            s += c->to_str() + ":";
        }
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
};

typedef std::shared_ptr<Structure> StructureOP;
typedef std::vector<StructureOP>   StructureOPs;

}

#endif //RNAMAKE_MOTIF_STATE_STRUCTURE_H
