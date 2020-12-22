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

enum ResType {
    ADE, CYT, GUA, URA, NONE
};

typedef std::vector<ResType> ResTypes;


ResType
convert_res_name_to_type(
        char);


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
    chain_id_(chain_id),
    dot_bracket_(dot_bracket),
    i_code_(i_code),
    name_(name),
    num_(num),
    uuid_(uuid) {
        res_type_ = convert_res_name_to_type(name_[0]);
    }
    
    inline
    Residue(
        Residue const & r):
    chain_id_(r.chain_id_),
    dot_bracket_(r.dot_bracket_),
    i_code_(r.i_code_),
    name_(r.name_),
    num_(r.num_),
    res_type_(r.res_type_),
    uuid_(r.uuid_)
    {}
    
    Residue(
        String const & s) {
        
        Strings spl = base::split_str_by_delimiter(s, ",");
        if(spl.size() < 4) {
            throw Exception("cannot build secondary_structure::Residue from str: " + s);
        }
        
        name_         = spl[0];
        dot_bracket_  = spl[1];
        num_          = std::stoi(spl[2]);
        chain_id_     = spl[3];
        uuid_         = util::Uuid();
        if(spl.size() == 5) {
            i_code_ = spl[4];
        }
        res_type_ = convert_res_name_to_type(name_[0]);
    }
    
    ~Residue() {}
    
public:

    inline
    String
    to_str() {
        std::stringstream ss;
        ss << name_ << "," << dot_bracket_ << "," << num_ << "," << chain_id_ << "," << i_code_;
        return ss.str();
    }
    
public: //getters
    
    inline
    String const &
    name() { return name_; }
    
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

    inline  // added by CJ for build_sqlite_libraries
    void
    i_code(String const& code) {
        i_code_= code;
    }

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
        name_ = name;
        res_type_ = convert_res_name_to_type(name_[0]);
    }
    

private:
    int num_;
    ResType res_type_;
    String name_, dot_bracket_, chain_id_, i_code_;
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
