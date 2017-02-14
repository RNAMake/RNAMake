//
// Created by Joseph Yesselman on 1/29/17.
//

#ifndef RNAMAKE_PRIMITIVES_BASEPAIR_H
#define RNAMAKE_PRIMITIVES_BASEPAIR_H

#include "util/uuid.h"
#include "primitives/residue.h"

namespace primitives {

class Basepair {
public:
    // WC: watson crick basepair
    // GU: gu basepair
    // NC : mismatched basepair (non-conical)
    enum BasepairType { WC, GU, NC};

public:
    inline
    Basepair(Uuid const & uuid):
            uuid_(uuid) {}

    virtual
    ~Basepair() {}

public:
    /**
    * equal operator checks whether the unique indentifier is the same
    * @param   other   another basepair to check if its the same
    */
    inline
    bool
    operator==(Basepair & other) {
        return uuid_ == other.uuid_;
    }

    inline
    bool
    operator!=(Basepair & other) const {
        return uuid_ != other.uuid_;
    }



protected:
    inline
    Basepair():
            uuid_(Uuid()) {}

public:
    inline
    Uuid const &
    uuid() { return uuid_; }

protected:
    Uuid uuid_;
};

template <typename Restype>
String
calc_bp_name(std::vector<std::shared_ptr<Restype>> const & res) {
    auto res1 = res[0];
    auto res2 = res[1];

    auto res1_name = String("");
    auto res2_name = String("");

    if(res1->i_code() == ' ') {
        res1_name = res1->chain_id()+std::to_string(res1->num());
    }
    else {
        res1_name = res1->chain_id()+std::to_string(res1->num())+res1->i_code();

    }

    if(res2->i_code() == ' ') {
        res2_name = res2->chain_id()+std::to_string(res2->num());
    }
    else {
        res2_name = res2->chain_id()+std::to_string(res2->num())+res2->i_code();
    }

    if(res1->chain_id() < res2->chain_id()) { return res1_name+"-"+res2_name; }
    if(res2->chain_id() < res1->chain_id()) { return res2_name+"-"+res1_name; }

    if(res1->num() < res2->num()) { return res1_name+"-"+res2_name; }
    else                          { return res2_name+"-"+res1_name; }

}

}

#endif //TEST_BASEPAIR_H
