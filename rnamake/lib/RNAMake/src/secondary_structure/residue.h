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

#include "base/types.h"
#include "base/string.h"
#include "util/uuid.h"

namespace sstruct {
    
class SecondaryStructureException : public std::runtime_error {
public:
    SecondaryStructureException(
        String const & message) :
    std::runtime_error(message)
    {}
    
};
    
    
class Residue {
public:
    inline
    Residue(
        String const & name,
        String const & dot_bracket,
        int const & num,
        String const & chain_id,
        Uuid const & uuid,
        String const & i_code=""):
    name_(name),
    dot_bracket_(dot_bracket),
    num_(num),
    chain_id_(chain_id),
    uuid_(uuid),
    i_code_(i_code)
    {}
    
    inline
    Residue(
        Residue const & r):
    name_(r.name_),
    dot_bracket_(r.dot_bracket_),
    num_(r.num_),
    chain_id_(r.chain_id_),
    uuid_(r.uuid_),
    i_code_(r.i_code_)
    {}
    
    Residue(
        String const & s) {
        
        Strings spl = split_str_by_delimiter(s, ",");
        name_         = spl[0];
        dot_bracket_  = spl[1];
        num_          = std::stoi(spl[2]);
        chain_id_     = spl[3];
        uuid_         = Uuid();
        if(spl.size() == 5) {
            i_code_ = spl[4];
        }

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
    
    inline
    Uuid const &
    uuid() { return uuid_; }

public: //setters
    
    inline
    void
    uuid(Uuid const & nuuid) { uuid_ = nuuid; }
    
    inline
    void
    name(String const & name) { name_ = name;} 
    

private:
    int num_;
    String name_, dot_bracket_, chain_id_, i_code_;
    Uuid uuid_;

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
    
} //sstruct

#endif /* defined(__RNAMake__ss_residue__) */
