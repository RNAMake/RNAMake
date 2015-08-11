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
    
    ~Residue() {}
    
public:
    
    inline
    Residue
    copy() {
        return Residue(name_, dot_bracket_, num_, chain_id_, uuid_, i_code_);
    }
    
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
    

private:
    int num_;
    String name_, dot_bracket_, chain_id_, i_code_;
    Uuid uuid_;

};
    
Residue
str_to_residue(String const & s);

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
