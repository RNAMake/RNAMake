//
//  secondary_structure_tree.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure_tree.h"
#include <queue>

int
get_brack_pair(
    int pos,
    String const & ss) {
    int bracket_count = 0;
    int i = -1;
    for(auto const & s : ss) {
        i++;
        if( i < pos) { continue; }
        if (s == '(') { bracket_count++; }
        if (s == ')') { bracket_count--; }
        if(bracket_count == 0) { return i; }
    }
    throw "could not find bracket pair";
    
}

int
get_dot_bounds(int, String const &, int) {
    return 1;
}


void
SecondaryStructureTree::_build_tree(
    String & ss,
    String & seq) {
    
    SecondaryStructureNode* node;
    if     ((ss[0] == '.' || ss.back() == '.') && nodes_.size() == 0) {
        
    }
    
    else if(ss[0] == '(') {
        node = new SSN_Basepair(ss, seq, NULL);
    }
    
    std::queue<SecondaryStructureNode*> open_nodes;
    SecondaryStructureNode* current;
    open_nodes.push(node);
    while (! open_nodes.empty()) {
        current = open_nodes.front();
        open_nodes.pop();
        for(auto const & c : current->children()) { open_nodes.push(c); }
        nodes_.push_back(current);
    }
    
    //std::cout << nodes_.size() << std::endl;
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SecondaryStructureNode
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SSandSeqOP
SecondaryStructureNode::get_ss_and_seq() {
    return SSandSeqOP(new SSandSeq("", ""));
}

void
SecondaryStructureNode::_assign_children(
    String & ss,
    String & seq) {
    
    if     (ss.back() == '.') {
        
    }
    
    else if(ss[0] == '(') {
        int pair = get_brack_pair(0, ss);
        String nss = ss.substr(0,pair+1);
        String nseq = seq.substr(0,pair+1);
        SecondaryStructureNode* child = new SSN_Basepair(nss, nseq, this);
        children_.push_back(child);
    }

}


SSN_Basepair::SSN_Basepair(
    String & ss,
    String & seq,
    SecondaryStructureNode* const & parent) {
    
    res1_ = seq.back();
    res2_ = seq[0];
    seq.pop_back();
    seq.erase(seq.begin(), seq.begin()+1);
    ss.pop_back();
    ss.erase(ss.begin(), ss.begin()+1);
    bp_type_ = String();
    bp_type_.push_back(res1_); bp_type_.push_back(res2_);
    _assign_children(ss, seq);

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
    String & ss,
    String & seq,
    SecondaryStructureNode* const & parent) {
    
}












