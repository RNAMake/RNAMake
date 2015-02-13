//
//  secondary_structure_tree.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__secondary_structure_tree__
#define __RNAMake__secondary_structure_tree__
#include <stdio.h>
#include <iostream>
#include "types.h"
#include "secondary_structure_tree.fwd.h"

struct SSandSeq {
public:
    SSandSeq(
        String const & nss,
        String const & nseq):
        ss ( nss ),
        seq ( nseq ) {}
    
public:
    String ss, seq;
};

typedef std::shared_ptr<SSandSeq> SSandSeqOP;

int
get_brack_pair(int, String const &);

int
get_dot_bounds(int, String const &, int);


class SecondaryStructureNode {
public:
    
    ~SecondaryStructureNode() {
        delete parent_;
    }
    
    virtual
    SSandSeqOP
    get_ss_and_seq();

public:
    
    inline
    SecondaryStructureNodes const &
    children() { return children_; }
    
protected:
    
    void
    _assign_children(
        String &,
        String &);
    
protected:
    SecondaryStructureNode* parent_;
    SecondaryStructureNodes children_;
    
};


class SSN_Basepair : public SecondaryStructureNode {
public:
    SSN_Basepair(
        String &,
        String &,
        SecondaryStructureNode* const &);
    
    SSandSeqOP
    get_ss_and_seq();

private:
    char res1_, res2_;
    String bp_type_, ss_type_;
    
    
};

class SSN_Bulge : public SecondaryStructureNode {
    SSN_Bulge(
        String &,
        String &,
        SecondaryStructureNode* const &);

};


class SecondaryStructureTree {
public:
    SecondaryStructureTree(
        String const & ss,
        String const & seq) {
        String css (ss);
        String cseq (seq);
        _build_tree(css, cseq);

    }
    
    ~SecondaryStructureTree() {
        for(int i = 0; i < nodes_.size(); i++) { delete nodes_[i]; }
    }
    
public:
    SSandSeqOP
    get_ss_and_seq() { return nodes_[0]->get_ss_and_seq(); }
    
private:
    
    void
    _build_tree(String &, String &);

private:
    SecondaryStructureNodes nodes_;
    
};

#endif /* defined(__RNAMake__secondary_structure_tree__) */
