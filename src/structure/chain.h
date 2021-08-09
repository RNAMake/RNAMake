//
// Created by Joseph Yesselman on 11/29/17.
//

#ifndef RNAMAKE_ALL_ATOM_CHAIN_H
#define RNAMAKE_ALL_ATOM_CHAIN_H

#include <primitives/chain.h>
#include <structure/residue.h>

namespace structure {

class Chain : public primitives::Chain<Residue> {
public:
    typedef primitives::Chain<Residue> ParentClass;

public:
    inline
    Chain(
            Residues const & residues): ParentClass(residues) {}

    inline
    Chain(
            Chain const & c) {
        residues_ = Residues();
        for(auto const & r : c) {
            residues_.push_back(Residue(r));
        }
    }

    Chain(
            String const & s,
            ResidueTypeSet const & rts) {

        residues_ = Residues();
        Strings spl = base::split_str_by_delimiter(s, ";");
        for(auto const & r_str : spl) {
            if (r_str.length() < 3) { continue; }
            residues_.push_back(Residue(r_str, rts));
        }
    }


    virtual
    ~Chain() {}

public:
    bool
    operator == (
            Chain const & c) const {
        return is_equal(c);
    }

    bool
    operator != (
            Chain const & c) const {
        return !(is_equal(c));
    }

public:
    bool
    is_equal(
            Chain const & c,
            bool check_uuid = true) const {

        if(residues_.size() != c.residues_.size()) { return false; }

        for(int i = 0; i < c.get_length(); i++) {
            if(!(residues_[i].is_equal(c.residues_[i], check_uuid))) { return false; }
        }
        return true;
    }


public: // non const methods

    void
    move(
            math::Point const & p) {
        for(auto & r : residues_) { r.move(p); }
    }

    void
    transform(
            math::Matrix const & r,
            math::Vector const & t,
            math::Point & dummy) {

        for(auto & res : residues_) { res.transform(r, t, dummy); }
    }

    inline
    void
    transform(
            math::Matrix const & r,
            math::Vector const & t) {
        auto dummy = math::Point();
        transform(r, t, dummy);
    }

public:

    inline
    String
    get_str() const {
        auto s = String("");
        for (auto const & r : residues_) { s += r.get_str() + ";"; }
        return s;
    }

    String
    get_pdb_str(
            int,
            int,
            char) const;

    inline
    String
    get_pdb_str(
            int acount = 1) const {
        return get_pdb_str(acount, get_first().get_num(), get_first().get_chain_id());
    }

    void
    write_pdb(
            String const & f_name) const {
        auto out = std::ofstream();
        out.open(f_name);
        out << get_pdb_str();
        out.close();
    }


};

typedef std::shared_ptr<Chain> ChainOP;
typedef std::vector<ChainOP>   ChainOPs;


}


#endif //RNAMAKE_ALL_ATOM_CHAIN_H
