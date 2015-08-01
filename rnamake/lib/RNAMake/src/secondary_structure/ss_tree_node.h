//
//  ss_tree_node.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__ss_tree_node__
#define __RNAMake__ss_tree_node__

#include <stdio.h>
#include "base/string.h"
#include "secondary_structure/ss_tree_node.fwd.h"
#include "secondary_structure/chain.h"

namespace sstruct {

class SS_NodeData {
public:
    enum SS_Type { SS_BP, SS_BULGE, SS_HAIRPIN, SS_NWAY, SS_PSEUDO_BP, SS_SEQ_BREAK, SS_MULTI };
    enum Bound_Side { LEFT, RIGHT };
    
    inline
    SS_NodeData(
        SS_Type type,
        ChainOPs const & ss_chains,
        ResidueOPs const & bounds):
    ss_chains_(ss_chains),
    type_(type),
    bounds_(bounds)
    {}
        
public: //getters
    
    ResidueOP const &
    bound_side(
        int pos,
        Bound_Side side) {
    
        try {
            auto ss_chain = ss_chains_.at(pos);
            if(ss_chain->length() == 0) {
                return bounds_[pos];
            }
            if(side == Bound_Side::LEFT) { return ss_chain->first(); }
            else {                        return ss_chain->last();  }
        }
        
        catch(std::out_of_range) {
            throw SecondaryStructureException("cannot access chain pos in bound_side");
        }
        catch(...) {
            throw std::runtime_error("unexpected error in bound_side");
        }
        
    }
    
    inline
    SS_Type const &
    type() { return type_; }
    
    inline
    String
    what() {
        if     (type_ == SS_BP)       { return "SS_BP";       }
        else if(type_ == SS_BULGE)    { return "SS_BULGE";    }
        else if(type_ == SS_HAIRPIN)  { return "SS_HAIRPIN";  }
        else if(type_ == SS_NWAY)     { return "SS_NWAY";     }
        else if(type_ == SS_PSEUDO_BP){ return "SS_PSEUDO_BP";}
        else if(type_ == SS_SEQ_BREAK){ return "SS_SEQ_BREAK";}
        else { throw std::runtime_error("unknown SS_TYPE"); }
    }
    
    inline
    String
    sequence() {
        String seq;
        int i = 0;
        for (auto const & c : ss_chains_) {
            seq += c->sequence();
            if(i+1 != ss_chains_.size()) { seq += "&"; }
            i++;
        }
        return seq;
    }
    
public: //getters
    
    inline
    ChainOPs
    ss_chains() { return ss_chains_; }
    
    inline
    ResidueOPs
    bounds() { return bounds_; }
    
    
protected:
    ChainOPs ss_chains_;
    ResidueOPs bounds_;
    SS_Type type_;
};



} //sstruct



#endif /* defined(__RNAMake__ss_tree_node__) */
