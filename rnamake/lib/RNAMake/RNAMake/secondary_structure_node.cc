//
//  secondary_structure_node.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure_node.h"
#include "FileIO.h"
#include <algorithm>
#include <iostream>


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SecondaryStructureNode
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SSandSeqOP
SecondaryStructureNode::get_ss_and_seq() {
    return SSandSeqOP(new SSandSeq("", ""));
}

int
SecondaryStructureNode::_assign_children(
    String const & ss,
    String const & seq,
    int xb,
    int yb) {
    
    if     (ss[yb-1] == '.' || ss[xb+1] == '.') {
        SecondaryStructureNodeOP child (new SSN_Bulge(ss, seq, xb+1, yb-1, shared_from_this()));
        children_.push_back(child);
        return yb;
    }
    
    else if(ss[xb+1] == '(') {
        int pair = get_brack_pair(xb+1, ss);
        SecondaryStructureNodeOP child (new SSN_Basepair(ss, seq, xb+1, pair, shared_from_this()));
        children_.push_back(child);
        if(yb - child->y_bounds() > 1) {
            return child->y_bounds();
        }
        return yb;
    }
    
    //finished
    return yb;
    
}

void
SecondaryStructureNode::assign_all_children(
    String const & ss,
    String const & seq) {
 
    //if hairpin has no children, continue
    if(ss_type_ == SSN_HAIRPIN) { return; }
    int xb = x_bounds();
    int yb = y_bounds();
    int current = xb;
    while (current != yb) {
        current = _assign_children(ss, seq, current, yb);
    }
}


SSN_Basepair::SSN_Basepair(
    String const & ss,
    String const & seq,
    int x_pos,
    int y_pos,
    SecondaryStructureNodeOP const & parent) {
    
    x_pos_ = x_pos;
    y_pos_ = y_pos;
    x_length_ = 1;
    y_length_ = 1;
    res1_ = seq[x_pos];
    res2_ = seq[y_pos];
    bp_type_ = String();
    bp_type_.push_back(res1_); bp_type_.push_back(res2_);
    ss_type_ = SSN_BP;
    children_ = SecondaryStructureNodeOPs();
    
}

SSandSeqOP
SSN_Basepair::get_ss_and_seq() {
    String ss, seq;
    for( auto const & c : children_) {
        SSandSeqOP ss_and_seq = c->get_ss_and_seq();
        ss += ss_and_seq->ss;
        seq += ss_and_seq->seq;
    }
    ss = "(" + ss + ")";
    seq = res1_ + seq + res2_;
    SSandSeqOP result ( new SSandSeq(ss, seq));
    return result;
}

SSN_Bulge::SSN_Bulge(
    String const & ss,
    String const & seq,
    int x_pos,
    int y_pos,
    SecondaryStructureNodeOP const & parent) {
    
    x_pos_ = x_pos;
    y_pos_ = y_pos;
    x_seq_ = String();
    y_seq_ = String();
    parent_ = parent;
    children_ = SecondaryStructureNodeOPs();
    ss_type_ = SSN_TWOWAY;
    x_length_ = 0;
    y_length_ = 0;
    int end = get_dot_bounds(x_pos, ss, 0);
    if(end-x_pos > -1 && ss[x_pos] == '.') {
        for(int i = x_pos; i <= end; i++) { x_seq_ += seq[i]; }
    }
    x_length_ = (int)x_seq_.length();
    //is hairpin
    if(end == y_pos_) { ss_type_ = SSN_HAIRPIN; return; }
    end = get_dot_bounds(y_pos, ss, 1);
    if(y_pos-end > -1 && ss[y_pos] == '.') {
        for(int i = end; i <= y_pos; i++) { y_seq_ += seq[i]; }
    }
    y_length_ = (int)y_seq_.length();
}

SSandSeqOP
SSN_Bulge::get_ss_and_seq() {
    String ss, seq;
    for( auto const & c : children_) {
        SSandSeqOP ss_and_seq = c->get_ss_and_seq();
        ss += ss_and_seq->ss;
        seq += ss_and_seq->seq;
    }
    
    for(int i = 0; i < x_seq_.length(); i++) { ss = '.' + ss; }
    for(int i = 0; i < y_seq_.length(); i++) { ss = ss + '.'; }
    seq = x_seq_ + seq + y_seq_;
    SSandSeqOP result ( new SSandSeq(ss, seq));
    return result;
}


SSN_Junction::SSN_Junction(
    SecondaryStructureNodeOP org) {
    ss_type_ = SSN_NWAY;
    parent_ = org->parent();
    
    //parent_->replace_child(org, shared_from_this());
    x_pos_ = org->x_pos();
    y_pos_ = org->y_pos();
    Strings end_seqs = split_str_by_delimiter(org->seq(), "+");
    seqs_ = Strings();
    seqs_.push_back(end_seqs[0]);
    for(auto const & c : org->children()) {
        if(c->ss_type() != SSN_TWOWAY) { children_.push_back(c); }
        else {
            seqs_.push_back(c->x_seq());
            if(c->children().size() > 1) { throw "too many children to make SSN_NWAY"; }
            for (auto const & c2 : c->children()) { children_.push_back(c2); }
        }
    }
    seqs_.push_back(end_seqs[1]);
}

SSandSeqOP
SSN_Junction::get_ss_and_seq() {
    
    String ss, seq;
    for(int i = 0; i < seqs_[0].length(); i++) { ss = '.' + ss; }
    seq = seqs_[0];
    int j = 1;
 
    for( auto const & c : children_) {
        SSandSeqOP ss_and_seq = c->get_ss_and_seq();
        ss += ss_and_seq->ss;
        seq += ss_and_seq->seq;
        
        for(int i = 0; i < seqs_[j].length(); i++) { ss += '.'; }
        seq += seqs_[j];
        //std::cout << seqs_[j] << std::endl;
        j++;
    }
    
    SSandSeqOP result ( new SSandSeq(ss, seq));
    return result;
}

