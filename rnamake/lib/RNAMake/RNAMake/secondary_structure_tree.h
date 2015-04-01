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
#include <algorithm>
#include <map>
#include "types.h"
#include "secondary_structure_node.fwd.h"
#include "secondary_structure_node.h"


class SecondaryStructureTree {
public:
    
    SecondaryStructureTree():
    nodes_( SecondaryStructureNodeOPs() )
    {}
    
    SecondaryStructureTree(
        String const & ss,
        String const & seq) {
        String css (ss);
        String cseq (seq);
        _build_tree(css, cseq);
        ss_ = ss;
        seq_ = seq;

        helices_ = std::vector<SecondaryStructureNodeOPs>();
        SecondaryStructureNodeOPs helix;
        for(auto const & n : nodes_) {
            if(n->ss_type() == SSN_BP) {
                helix.push_back(n);
            }
            else {
                helices_.push_back(helix);
                helix = SecondaryStructureNodeOPs();
            }
        }
        if(helix.size() > 0) { helices_.push_back(helix); }
    }
    
    ~SecondaryStructureTree() {}
    
public:
    SSandSeqOP
    get_ss_and_seq() const { return nodes_[0]->get_ss_and_seq(); }
   
    inline
    SecondaryStructureNodeOPs
    get_bulges() const {
        SecondaryStructureNodeOPs bulges;
        for (auto const & n : nodes_) {
            if(n->ss_type() == SSN_TWOWAY) { bulges.push_back(n); }
        }
        return bulges;
    }
    
    inline
    SecondaryStructureNodeOPs
    get_designable_bps() const {
        SecondaryStructureNodeOPs bps;
        
        for ( auto const & n : nodes_) {
            if(n->ss_type() == SSN_BP && n->bp_type().compare("NN") == 0) {
                bps.push_back(n);
            }
        }
        
        return bps;
    }
    
    int
    gc_pairs() const {
        int count = 0;
        for (auto const & n : nodes_ ) {
            if(n->ss_type() == SSN_BP && (n->bp_type_num() == 0 || n->bp_type_num() == 1)) {
                count++;
            }
        }
        
        return count;
    }
    
    int
    gu_pairs() const {
        int count = 0;
        for (auto const & n : nodes_ ) {
            if(n->ss_type() == SSN_BP && (n->bp_type_num() == 4 || n->bp_type_num() == 5)) {
                count++;
            }
        }
        
        return count;
    }
    
    int
    ua_pairs() const {
        int count = 0;
        for (auto const & n : nodes_ ) {
            if(n->ss_type() == SSN_BP && (n->bp_type_num() == 2 || n->bp_type_num() == 3)) {
                count++;
            }
        }
        
        return count;
    }
    
    std::vector<SecondaryStructureNodeOPs> const &
    helices() const { return helices_; }
    
    std::map<int, int>
    get_pairmap() const {
        std::map<int, int> pairmap;
        
        for(auto const & n : nodes_) {
            if(n->ss_type() != SSN_BP) { continue; }
            pairmap[n->x_pos()+1] = n->y_pos()+1;
            pairmap[n->y_pos()+1] = n->x_pos()+1;
        }
        
        return pairmap;
    }
    
    inline
    String const &
    ss() const {
        return ss_;
    }
    
    String const
    seq() const {
        String seq = seq_;
        for(auto const & n : nodes_) {
            if(n->ss_type() != SSN_BP) { continue;}
            seq[ n->x_pos() ] = n->res1();
            seq[ n->y_pos() ] = n->res2();
        }
        return seq;
    }
    
    
private:
    
    void
    _build_tree(String &, String &);

private:
    SecondaryStructureNodeOPs nodes_;
    std::vector<SecondaryStructureNodeOPs> helices_;
    String ss_, seq_;
    
};

#endif /* defined(__RNAMake__secondary_structure_tree__) */
