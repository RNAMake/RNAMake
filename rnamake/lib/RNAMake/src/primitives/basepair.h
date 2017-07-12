//
// Created by Joseph Yesselman on 1/29/17.
//

#ifndef RNAMAKE_PRIMITIVES_BASEPAIR_H
#define RNAMAKE_PRIMITIVES_BASEPAIR_H

#include "util/uuid.h"
#include "util/x3dna.h"
#include "primitives/residue.h"

/*
 * Exception for basepair
 */
class BasepairException : public std::runtime_error {
public:
    /**
     * Standard constructor for BasepairException
     * @param   message   Error message for basepair
     */
    BasepairException(String const & message) :
            std::runtime_error(message) {}
};

namespace primitives {

class Basepair {
public:
    // WC: watson crick basepair
    // GU: gu basepair
    // NC : mismatched basepair (non-conical)
    enum BasepairType { WC, GU, NC};

public:
    inline
    Basepair(
            Uuid const & res1_uuid,
            Uuid const & res2_uuid,
            Uuid const & uuid):
            res1_uuid_(res1_uuid),
            res2_uuid_(res2_uuid),
            uuid_(uuid) {}

    virtual
    ~Basepair() {}

public:
    Uuid const &
    partner(Uuid const &);

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
            res1_uuid_(Uuid()),
            res2_uuid_(Uuid()),
            uuid_(Uuid()) {}

public:
    inline
    Uuid const &
    uuid() { return uuid_; }

    inline
    Uuid const &
    res1_uuid() const { return res1_uuid_; }

    inline
    Uuid const &
    res2_uuid() const { return res2_uuid_; }


protected:
    Uuid uuid_;
    Uuid res1_uuid_, res2_uuid_;
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

template <typename Restype>
Basepair::BasepairType
get_bp_type(std::vector<std::shared_ptr<Restype>> const & res,
            X3dna::X3dnaBPType const & x3dna_bp_type) {
    auto bp_str = String();
    bp_str += res[0]->name(); bp_str += res[1]->name();
    auto wc_names = Strings{"GC", "CG", "AU", "UA"};
    // not in WC orientation, is non-cononical
    if(x3dna_bp_type != X3dna::X3dnaBPType::cWUW) { return Basepair::BasepairType::NC; }
    // has correct residues to be a WC basepair
    if (std::find(wc_names.begin(), wc_names.end(), bp_str) != wc_names.end()) { return Basepair::BasepairType::WC; }
    if (bp_str == "GU" || bp_str == "UG" ) { return Basepair::BasepairType::GU; }
    return Basepair::BasepairType::NC;
}

}


#endif //TEST_BASEPAIR_H
