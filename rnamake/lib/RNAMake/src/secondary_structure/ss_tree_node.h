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


class SS_NodeData {
public:
    enum SS_Type { SS_BP, SS_BULGE, SS_HAIRPIN, SS_NWAY, SS_PSEUDO_BP, SS_SEQ_BREAK };
    enum Bound_Side { LEFT, RIGHT };
    
    inline
    SS_NodeData(
        Strings const & seq,
        SS_Type type,
        std::vector<Ints> bounds):
    seq_(seq),
    type_(type),
    bounds_(bounds)
    {}
        
public: //getters
    
    Ints const &
    bounds(
        int pos) {
        if(pos >= bounds_.size()) {
            throw std::runtime_error("cannot call bounds in SS_TreeNode, with this pos");
        }
        
        return bounds_[pos];
    }

    int const &
    bound_side(
        int pos,
        Bound_Side side) {
        
        Ints bound = bounds(pos);
        return bound[(int)side];
    }
    
    inline
    SS_Type const &
    type() { return type_; }
    
    inline
    int
    nstrands() { return (int)seq_.size(); }
    
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
    Strings
    seqs() { return seq_; }
    
public: //virtual getters
    
    virtual
    inline
    String
    seq() {
        String seq;
        int i =-1;
        for(auto const & s : seq_) {
            i++;
            seq += s;
            if(i < seq_.size()-1) { seq += "+"; }
        }
        return seq;
    }
    
protected:
    Strings seq_;
    std::vector<Ints> bounds_;
    SS_Type type_;
};

class SS_NodeDataBP : public SS_NodeData {
public:
    SS_NodeDataBP(
        Strings const & s,
        SS_Type ntype,
        std::vector<Ints> bounds):
    SS_NodeData(s, ntype, bounds)
    { }
    
    
public:
    
    inline
    String
    seq() { return seq_[0]+seq_[1]; }
    
};




#endif /* defined(__RNAMake__ss_tree_node__) */
