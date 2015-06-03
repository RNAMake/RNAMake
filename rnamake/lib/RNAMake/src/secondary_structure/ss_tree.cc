//
//  ss_tree.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <queue>
#include <iostream>

#include "secondary_structure/ss_tree.h"
#include "secondary_structure/ss_tree_node.h"

SS_Tree::SS_Tree(
    String const & ss,
    String const & seq):
    ss_(ss),
    seq_(seq),
    tree_(Tree<SS_NodeDataOP>()) {
        
    _build_tree();
}

void
SS_Tree::_build_tree() {
    int xb = 0, yb = (int)seq_.size()-1;
    SS_NodeDataOP current = _assign_new_node(xb, yb);
    int index;
    
    using NodeandIndex = std::pair<SS_NodeDataOP, int>;
    std::queue<NodeandIndex> open_nodes;
    open_nodes.push(NodeandIndex(current, -1));
    
    while(! open_nodes.empty()) {
        NodeandIndex current_pair = open_nodes.front();
        current = current_pair.first;
        
        if(current->type() == SS_NodeData::SS_Type::SS_HAIRPIN) {
            index = tree_.add_data(current, 0, current_pair.second);
            open_nodes.pop();
            continue;
        }
        
        int xb = current->bound_side(0, SS_NodeData::Bound_Side::RIGHT)+1;
        int yb = current->bound_side(1, SS_NodeData::Bound_Side::LEFT) -1;
        
        std::vector<SS_NodeDataOP> next_level = _build_tree_level(xb, yb);
        
        // Check for multi way junctions will be defined as a bulge as a child of another bulge
        if(current->type() == SS_NodeData::SS_Type::SS_BULGE) {
            std::vector<SS_NodeDataOP> part_of_nway;
            std::vector<SS_NodeDataOP> not_part_of_nway;
            for(auto const & n : next_level) {
                if(n->type() == SS_NodeData::SS_Type::SS_BULGE ||
                   n->type() == SS_NodeData::SS_Type::SS_HAIRPIN) { part_of_nway.push_back(n); }
                else { not_part_of_nway.push_back(n); }
            }
            next_level = not_part_of_nway;
            
            if(part_of_nway.size() > 0) {
                Strings seq;
                std::vector<Ints> bounds;
                
                seq.push_back(current->seqs()[0]);
                bounds.push_back(current->bounds(0));

                seq.push_back(current->seqs()[1]);
                bounds.push_back(current->bounds(1));
                
                for(auto const & n : part_of_nway) {
                    Strings child_seq = n->seqs();
                    for(int i = 0; i < child_seq.size(); i++) {
                        if(child_seq[i].size() > 0) {
                            seq.push_back(child_seq[i]);
                            bounds.push_back(n->bounds(i));
                        }
                    }
                    
                    if(n->type() == SS_NodeData::SS_Type::SS_BULGE) {
                        int nxb = n->bound_side(0, SS_NodeData::Bound_Side::RIGHT)+1;
                        int nyb = n->bound_side(1, SS_NodeData::Bound_Side::LEFT) -1;
                        std::vector<SS_NodeDataOP> next_level_2 = _build_tree_level(nxb, nyb);
                        for(auto const & n2 : next_level_2) { next_level.push_back(n2); }
                    }
                }
                
                
                current = std::make_shared<SS_NodeData>(seq, SS_NodeData::SS_Type::SS_NWAY, bounds);
            }
            
        }
        
        index = tree_.add_data(current, 0, current_pair.second);
        open_nodes.pop();
        
        for(auto const & n : next_level) {
            open_nodes.push(NodeandIndex(n, index));
        }
    }
    
}

std::vector<SS_NodeDataOP>
SS_Tree::_build_tree_level(
    int xb,
    int yb) {
    
    std::vector<SS_NodeDataOP> next_level;
    while(xb <= yb) {
        SS_NodeDataOP child = _assign_new_node(xb, yb);
        xb = child->bound_side(1, SS_NodeData::Bound_Side::RIGHT)+1;
        next_level.push_back(child);
    }
    
    return next_level;
    
}


SS_NodeDataOP
SS_Tree::_assign_new_node(
    int xb,
    int yb) {
    
    Strings seq(2);
    std::vector<Ints> bounds(2);
    bounds[0] = Ints({xb-1, xb-1});
    bounds[1] = Ints({yb+1, yb+1});
    int hairpin = 0;
    SS_NodeDataOP current;
    if(ss_[xb] == '.') {
        int end_x = _get_dot_bounds(xb, 0);
        for(int i = xb; i <= end_x; i++) { seq[0] += seq_[i]; }
        if(end_x == yb) { hairpin = 1; }
        bounds[0] = Ints({xb, end_x});
    }

    if(ss_[yb] == '.' && hairpin == 0) {
        int end_y = _get_dot_bounds(yb, 1);
        for(int i = end_y; i <= yb; i++) { seq[1] += seq_[i]; }
        bounds[1] = Ints({end_y, yb});
    }
    
    if(ss_[xb] == '&' || ss_[xb] == '+') {
        bounds[0] = Ints({xb, xb});
        current = std::make_shared<SS_NodeData>(seq, SS_NodeData::SS_Type::SS_SEQ_BREAK, bounds);
        return current;
    }
    
    if(seq[0].length() == 0 && seq[1].length() == 0) {
        int pair = _get_brack_pair(xb);
        seq[0] = seq_[xb];
        seq[1] = seq_[pair];
        bounds[0] = Ints({xb, xb});
        bounds[1] = Ints({pair, pair});

        
        if(ss_[xb] == '(') {
            current = std::make_shared<SS_NodeDataBP>(seq, SS_NodeData::SS_Type::SS_BP, bounds);
        }
        else               {
            current = std::make_shared<SS_NodeDataBP>(seq, SS_NodeData::SS_Type::SS_PSEUDO_BP, bounds);
        }

    }
    else if(hairpin) {
        current = std::make_shared<SS_NodeData>(seq, SS_NodeData::SS_Type::SS_HAIRPIN, bounds);
    }

    else {
        current = std::make_shared<SS_NodeData>(seq, SS_NodeData::SS_Type::SS_BULGE, bounds);
    }
    
    return current;

}


SeqSS
SS_Tree::seq_from_nodes(
    SS_Nodes const & nodes) const {
    
    std::map<int, int> seen;
    for(auto const & n : nodes) {
        for(auto const & bound : n->data()->bounds()) {
            for(int i = bound[0]; i <= bound[1]; i++) {
                seen[i] = 1;
            }
        }
    }
    
    int start = -1;
    int range = 0;
    int found = 0;
    Strings seqs, sss;
    for(int i = 0; i < seq_.length(); i++) {
        if(seen.find(i) != seen.end()) {
            if(found == 0) {
                found = 1;
                start = i;
                range = 1;
            }
            else {
                range += 1;
            }
        }
        
        else {
            found = 0;
            String seq = seq_.substr(start, range);
            String ss  = ss_.substr(start, range);
            seqs.push_back(seq);
            sss.push_back(ss);
            
        }
    }
    
    if(found == 1) {
        String seq = seq_.substr(start, range);
        String ss  = ss_.substr(start, range);
        seqs.push_back(seq);
        sss.push_back(ss);
    }
    
    return SeqSS(seqs, sss);
    
}

int
SS_Tree::_get_brack_pair(
    int pos) {
    
    int bracket_count = 0;
    int i = -1;
    for(auto const & s : ss_) {
        i++;
        if( i < pos) { continue; }
        if (s == '(' || s == '[') { bracket_count++; }
        if (s == ')' || s == ']') { bracket_count--; }
        if(bracket_count == 0) { return i; }
    }
    
    throw std::runtime_error("could not find bracket pair");

}

int
SS_Tree::_get_dot_bounds(
    int pos,
    int reverse) {

    if(!reverse) {
        for(int i = pos+1; i < ss_.length(); i++) {
            if(ss_[i] != '.') { return i-1; }
        }
        return (int)ss_.length()-1;
    }
    else {
        for(int i = pos-1; i >= 0; i--) {
            if(ss_[i] != '.') { return i+1; }
        }
        return 0;
    }
    return pos;
}


