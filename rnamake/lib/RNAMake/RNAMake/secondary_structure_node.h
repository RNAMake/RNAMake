//
//  secondary_structure_node.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__secondary_structure_node__
#define __RNAMake__secondary_structure_node__

#include <stdio.h>
#include "secondary_structure_node.fwd.h"
#include "types.h"

enum SSN_Type {
    SSN_BP      = 0,
    SSN_TWOWAY  = 1,
    SSN_NWAY    = 2,
    SSN_HAIRPIN = 3
};


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
        //delete parent_;
    }
    
    virtual
    SSandSeqOP
    get_ss_and_seq();
    
public:
    
    void
    assign_all_children(
        String const &,
        String const &);
    
    inline
    virtual
    String
    seq() { return ""; }
    
    void
    replace_child(SecondaryStructureNode* & org_child,
                  SecondaryStructureNode* new_child) {
        int pos = (int)(std::find(children_.begin(), children_.end(), org_child) - children_.begin());
        children_[pos] = new_child;
    }
    
    

    
public:
    
    inline
    SecondaryStructureNodes const &
    children() { return children_; }
    
    inline
    int
    x_bounds() { return x_pos_ + x_length_ - 1; }
    
    inline
    int
    y_bounds() { return y_pos_ - y_length_ + 1; }
    
    inline
    SSN_Type const &
    ss_type() { return ss_type_; }
    
    inline
    SecondaryStructureNode*
    parent() { return parent_; }
    
    inline
    int
    x_pos() { return x_pos_; }
    
    inline
    int
    y_pos() { return y_pos_; }
    
    virtual
    inline
    String
    x_seq() { return ""; }
    
    
protected:
    
    int
    _assign_children(
        String const &,
        String const &,
        int,
        int);
    
protected:
    SecondaryStructureNode* parent_;
    SecondaryStructureNodes children_;
    SSN_Type ss_type_;
    int x_pos_, y_pos_, x_length_, y_length_;
    
    
};

class SSN_Basepair : public SecondaryStructureNode {
public:
    SSN_Basepair(
        String const &,
        String const &,
        int x_pos,
        int y_pos,
        SecondaryStructureNode* const &);
    
    SSandSeqOP
    get_ss_and_seq();
    
private:
    char res1_, res2_;
    String bp_type_;
    
    
};

class SSN_Bulge : public SecondaryStructureNode {
public:
    SSN_Bulge(
        String const &,
        String const &,
        int x_pos,
        int y_pos,
        SecondaryStructureNode* const &);
    
    SSandSeqOP
    get_ss_and_seq();

public:
    inline
    String
    seq() { return x_seq_ + "+" + y_seq_; }
    
    inline
    String
    x_seq() { return x_seq_; }
    
    
private:
    String x_seq_, y_seq_, bulge_type_;
};

class SSN_Junction : public SecondaryStructureNode {
public:
    SSN_Junction(
        SecondaryStructureNode*);
    
public:
    
    SSandSeqOP
    get_ss_and_seq();
    
private:
    Strings seqs_;
};

#endif /* defined(__RNAMake__secondary_structure_node__) */
