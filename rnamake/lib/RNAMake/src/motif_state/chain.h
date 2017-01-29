//
// Created by Joseph Yesselman on 1/28/17.
//

#ifndef RNAMAKE_MOTIF_STATE_CHAIN_H
#define RNAMAKE_MOTIF_STATE_CHAIN_H

#include "primitives/chain.h"
#include "motif_state/residue.h"

namespace state {

class Chain : public primitives::Chain<Residue>{
public:
    Chain(
            ResidueOPs const & residues):
            primitives::Chain<Residue>(residues) {}

    Chain(
            Chain const & c):
            primitives::Chain<Residue>() {

        residues_ = ResidueOPs(c.residues_.size());
        int i = 0;
        for (auto const & r : c.residues_) {
            residues_[i] = std::make_shared<Residue>(*r);
            i++;
        }
    }

    Chain(
            String const & s):
            primitives::Chain<Residue>() {

        residues_ = ResidueOPs();
        Strings spl = split_str_by_delimiter(s, ";");
        for (auto const & r_str : spl) {
            auto r = std::make_shared<Residue>(r_str);
            residues_.push_back(r);
        }
    }

    ~Chain() {}

public:
    String
    to_str();

public:

    inline
    void
    move(Point const & p) {
        for(auto & r : residues_) { r->move(p); }
    }

    inline
    void
    transform(Transform const & t) {
        for(auto & r : residues_) { r->transform(t); }
    }

    inline
    void
    fast_transform(
            Matrix const & r,
            Vector const & t) {
        for(auto & res : residues_) { res->fast_transform(r, t); }
    }

};

typedef std::shared_ptr<Chain> ChainOP;
typedef std::vector<ChainOP> ChainOPs;

}
#endif //RNAMAKE_MOTIF_STATE_CHAIN_H
