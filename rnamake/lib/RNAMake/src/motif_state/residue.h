//
// Created by Joseph Yesselman on 1/28/17.
//

#ifndef RNAMAKE_MOTIF_STATE_RESIDUE_H
#define RNAMAKE_MOTIF_STATE_RESIDUE_H

#include <sstream>

#include "primitives/residue.h"
#include "util/bead.h"

namespace state {

class Residue : public primitives::Residue {
public:
    inline
    Residue(
            char name,
            int num,
            char chain_id,
            char i_code,
            Beads const & beads):
            primitives::Residue(name, num, chain_id, i_code),
            beads_(beads) {}

    inline
    Residue(
            char name,
            int num,
            char chain_id,
            char i_code,
            Beads const & beads,
            Uuid const & uuid):
            primitives::Residue(name, num, chain_id, i_code, uuid),
            beads_(beads) {}

    inline
    Residue(
            Residue const & r,
            int new_uuid = 0):
            primitives::Residue(r.name_, r.num_, r.chain_id_, r.i_code_, r.uuid_),
            beads_(r.beads_) {

        if(new_uuid) { uuid_ = Uuid(); }
    }

    inline
    Residue(
            String const & s):
            primitives::Residue() {

        auto spl = split_str_by_delimiter(s, ",");
        name_     = spl[0][0];
        num_      = std::stoi(spl[1]);
        chain_id_ = spl[2][0];
        i_code_   = spl[3][0];
        uuid_     = Uuid();
        beads_    = Beads();
        for(int i = 4; i < spl.size()-1; i+=2) {
            auto b = Bead(spl[i] + "," + spl[i+1]);
            beads_.push_back(b);
        }

    }

public:
    String
    to_str() {
        auto ss = std::stringstream();
        ss << name_ << "," << std::to_string(num_) << "," << chain_id_ << "," << i_code_ << ",";
        for(auto const & b : beads_) {
            ss << b.to_str() << ",";
        }

        return ss.str();
    }

public:
    inline
    void
    move(Point const & p) {
        for(auto & b : beads_) { b.move(p); }
    }

    inline
    void
    transform(Transform const & t) {
        for(auto & b : beads_) { b.transform(t); }
    }

    inline
    void
    fast_transform(
            Matrix const & r,
            Vector const & t,
            Point & dummy) {
        for(auto & b : beads_) { b.fast_transform(r, t, dummy); }
    }

public: //getters

    inline
    Beads const &
    beads() { return beads_; }

    size_t
    num_beads() { return beads_.size(); }


private:
    Beads beads_;
};

typedef std::shared_ptr<Residue> ResidueOP;
typedef std::vector<ResidueOP>   ResidueOPs;

}
#endif //RNAMAKE_MOTIF_STATE_RESIDUE_H
