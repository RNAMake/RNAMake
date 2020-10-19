//
//  residue.h
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__ss_residue__
#define __RNAMake__ss_residue__

#include <stdio.h>
#include <sstream>
#include <memory>
#include <stdexcept>
#include <cassert>

#include "base/types.h"
#include "base/string.h"
#include "util/uuid.h"

namespace secondary_structure {

class Exception : public std::runtime_error {
public:
    Exception(
            String const & message) :
            std::runtime_error(message) {}
};

//https://en.wikipedia.org/wiki/Nucleic_acid_notation
enum class ResType {
    A, // A
    C, // C
    G, // G
    U, // U
    W, // "Weak" A or U
    S, // "Strong" C or G
    M, // "aMino" A or C
    K, // "Keto" G or U
    R, // "puRine", A or G
    Y, // "pYrimidine, C or U
    B, // not A, C G or U
    D, // not C, A G or U
    H, // not G, A C or G
    V, // not U, A C or G
    N  // A, C, G or U
};

typedef std::vector<ResType> ResTypes;


ResType
convert_res_name_to_type(
        char);

String
convert_res_type_to_str(
        ResType);

bool
does_restype_satisfy_constraint(
        ResType,
        ResType);

bool
is_restype_a_weak(
        ResType);

bool
is_restype_a_strong(
        ResType);

bool
is_restype_a_amino(
        ResType);

bool
is_restype_a_keto(
        ResType);

bool
is_restype_a_purine(
        ResType);

bool
is_restype_a_pyrimidine(
        ResType);

bool
is_restype_not_A(
        ResType);

bool
is_restype_not_C(
        ResType);

bool
is_restype_not_G(
        ResType);

bool
is_restype_not_U(
        ResType);

bool
is_restype_a_ambiguous_code(
        ResType);


class Residue {
public:
    inline
    Residue(
        String const & name,
        String const & dot_bracket,
        int const & num,
        String const & chain_id,
        util::Uuid const & uuid,
        String const & i_code=""):
    dot_bracket_(dot_bracket),
    num_(num),
    chain_id_(chain_id),
    uuid_(uuid),
    i_code_(i_code),
    res_type_(convert_res_name_to_type(name[0]))
    {}
    
    inline
    Residue(
        Residue const & r):
    dot_bracket_(r.dot_bracket_),
    num_(r.num_),
    chain_id_(r.chain_id_),
    uuid_(r.uuid_),
    i_code_(r.i_code_),
    res_type_(r.res_type_)
    {}
    
    Residue(
        String const & s) {
        
        Strings spl = base::split_str_by_delimiter(s, ",");
        if(spl.size() < 4) {
            throw Exception("cannot build secondary_structure::Residue from str: " + s);
        }

        dot_bracket_  = spl[1];
        num_          = std::stoi(spl[2]);
        chain_id_     = spl[3];
        uuid_         = util::Uuid();
        if(spl.size() == 5) {
            i_code_ = spl[4];
        }
        res_type_ = convert_res_name_to_type(spl[0][0]);
    }
    
    ~Residue() {}
    
public:

    inline
    String
    to_str() {
        std::stringstream ss;
        ss << name() << "," << dot_bracket_ << "," << num_ << "," << chain_id_ << "," << i_code_;
        return ss.str();
    }
    
public: //getters
    
    inline
    String
    name() { return convert_res_type_to_str(res_type_); }
    
    inline
    String const &
    dot_bracket() { return dot_bracket_; }
    
    inline
    int const &
    num() { return num_; }
    
    inline
    String const &
    chain_id() { return chain_id_; }
    
    inline
    String const &
    i_code() { return i_code_; }
    
    inline
    util::Uuid const &
    uuid() { return uuid_; }

    inline
    ResType
    res_type() { return res_type_; }
    
public: //setters
    
    inline
    void
    uuid(util::Uuid const & nuuid) { uuid_ = nuuid; }
    
    inline
    void
    name(String const & name) {
        res_type_ = convert_res_name_to_type(name[0]);
    }

    inline
    void
    res_type(ResType type) {
        res_type_ = type;
    }
    

private:
    int num_;
    ResType res_type_;
    String dot_bracket_, chain_id_, i_code_;
    util::Uuid uuid_;

};
    
typedef std::shared_ptr<Residue> ResidueOP;
typedef std::vector<ResidueOP> ResidueOPs;
    
struct res_less_than_key {
    inline
    bool
    operator() (ResidueOP const & r1, ResidueOP const & r2) {
        return (r1->num() < r2->num());
    }
};
    
} //secondary_structure

#endif /* defined(__RNAMake__ss_residue__) */
