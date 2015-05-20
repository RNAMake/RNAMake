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

struct SS_TreeNodeProto {
    SS_TreeNodeProto(
        SS_TreeNodeProtoOP const & nparent,
        int const nxb,
        int const nyb):
    
    parent(nparent),
    xb(nxb),
    yb(nyb),
    children(SS_TreeNodeProtoOPs())
    {}

    
    SS_TreeNodeProtoOP parent;
    SS_TreeNodeProtoOPs children;
    int xb, yb;
    
};


class SS_TreeNode {
public: //getters
    
    inline
    int
    xb() { return xb_; }
    
    inline
    int
    yb() { return yb_; }
    
    
public: //only bp getters
    inline
    virtual
    String
    bp_type() { throw "cannot call bp_type"; }
    
    inline
    virtual
    int
    bp_type_num() { throw "cannot call bp_type_num"; }
    
protected:
    Strings seq_;
    int xb_, yb_;
};


class SS_TreeNodeBasepair : public SS_TreeNode {
public:
    SS_TreeNodeBasepair(
        Strings const & s,
        int xb,
        int yb) {
        
        seq_ = s;
        xb_ = xb;
        yb_ = yb;
    }
    
    
public:
    inline
    virtual
    String
    bp_type() { return bp_type_; }
    
private:
    String bp_type_;
};




#endif /* defined(__RNAMake__ss_tree_node__) */
