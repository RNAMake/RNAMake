//
//  basepair.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__ss_basepair__
#define __RNAMake__ss_basepair__

#include <stdio.h>
#include <sstream>

//RNAMake Headers
#include "secondary_structure/residue.h"

namespace secondary_structure {

class Basepair {
public:
    inline
    Basepair(
        ResidueOP const & res1,
        ResidueOP const & res2,
        util::Uuid const & uuid):
    res1_(res1),
    res2_(res2),
    uuid_(uuid)
    {}
    
    ~Basepair() {}
    
public:
    
    inline
    String
    name() {
        std::stringstream ss;
        ss << res1_->chain_id() << res1_->num() << res1_->i_code();
        String str1 = ss.str();
        ss.str("");
        ss << res2_->chain_id() << res2_->num() << res2_->i_code();
        String str2 = ss.str();
        if(str1 < str2) { return str1+"-"+str2; }
        else            { return str2+"-"+str1; }
        
    }
    
    inline
    ResidueOP
    partner(
        ResidueOP const & r) {
        
        if     (r == res1_) { return res2_; }
        else if(r == res2_) { return res1_; }
        else {
            throw std::runtime_error("called partner with a resiude not in basepair in secondary_structure");
        }
        
    }
    
public:
    
    inline
    ResidueOP &
    res1() { return res1_; }
    
    inline
    ResidueOP &
    res2() { return res2_; }
    
    inline
    ResidueOP const &
    res1() const { return res1_; }
    
    inline
    ResidueOP const &
    res2() const { return res2_; }
    
    inline
    util::Uuid const &
    uuid() { return uuid_; }
    
private:
    ResidueOP res1_, res2_;
    util::Uuid uuid_;
    
};

typedef std::shared_ptr<Basepair> BasepairOP;
typedef std::vector<BasepairOP>   BasepairOPs;

inline
bool
is_gc_pair(BasepairOP const & bp) {
    if     (bp->res1()->res_type() == ResType::GUA && bp->res2()->res_type() == ResType::CYT) { return true; }
    else if(bp->res2()->res_type() == ResType::CYT && bp->res1()->res_type() == ResType::GUA) { return true; }
    else { return false; }
}
    
    
inline
bool
is_au_pair(BasepairOP const & bp) {
    if     (bp->res1()->res_type() == ResType::ADE && bp->res2()->res_type() == ResType::URA) { return true; }
    else if(bp->res2()->res_type() == ResType::URA && bp->res1()->res_type() == ResType::ADE) { return true; }
    else { return false; }
}
    
inline
bool
is_gu_pair(BasepairOP const & bp) {
    if     (bp->res1()->res_type() == ResType::GUA && bp->res2()->res_type() == ResType::URA) { return true; }
    else if(bp->res2()->res_type() == ResType::URA && bp->res1()->res_type() == ResType::GUA) { return true; }
    else { return false; }
}
    
    
}

#endif /* defined(__RNAMake__basepair__) */
