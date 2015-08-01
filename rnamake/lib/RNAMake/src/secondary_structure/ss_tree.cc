//
//  ss_tree.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <queue>
#include <iostream>

#include "secondary_structure/residue.h"
#include "secondary_structure/ss_tree.h"
#include "secondary_structure/ss_tree_node.h"

namespace sstruct {

struct NodeandIndex {
    SS_NodeDataOP node;
    int index;
    NodeandIndex(SS_NodeDataOP const & n, int i):
        index(i), node(n) {}
};
    
SS_Tree::SS_Tree(
    String const & seq,
    String const & ss):
    tree_(TreeDynamic<SS_NodeDataOP>())
{
    
    ss_ = SecondaryStructure(seq, ss);
    auto pad_residue1 = std::make_shared<Residue>("N", ".", 0, "N", Uuid());
    residues_ = ResidueOPs();
    residues_.push_back(pad_residue1);
    for(auto const & r : ss_.residues()) { residues_.push_back(r); }
    auto pad_residue2 = std::make_shared<Residue>("N", ".", ss_.residues().size(), "N", Uuid());
    residues_.push_back(pad_residue2);
    
    _build_tree();
}

void
SS_Tree::_build_tree() {
    int xb = 1, yb = (int)residues_.size()-2;
    auto ss_chains = std::vector<ChainOP>{ std::make_shared<Chain>(),
                                           std::make_shared<Chain>() };
    auto bounds = ResidueOPs { residues_[xb-1], residues_[yb+1] };
    auto current = std::make_shared<SS_NodeData>( SS_NodeData::SS_Type::SS_SEQ_BREAK,
                                                  ss_chains, bounds);

    std::queue<NodeandIndex> open_nodes;
    open_nodes.push(NodeandIndex(current, -1));
    int index = 0;
    while(! open_nodes.empty()) {
        NodeandIndex current_pair = open_nodes.front();
        open_nodes.pop();
        
        if(current_pair.node->type() == SS_NodeData::SS_Type::SS_HAIRPIN) {
            tree_.add_data(current_pair.node, current_pair.index);
            continue;
        }
        
        xb = _map_back_to_index(current_pair.node->bound_side(0, SS_NodeData::Bound_Side::RIGHT))+1;
        yb = _map_back_to_index(current_pair.node->bound_side(1, SS_NodeData::Bound_Side::LEFT ))-1;
       
        if(xb > yb) {
            index = tree_.add_data(current_pair.node, current_pair.index);
            auto seq_break = _check_from_chain_ends(xb-1, yb+1);
            if(seq_break != nullptr) { tree_.add_data(seq_break, index); }
            continue;
        }
        
        auto next_level = _build_tree_level(xb, yb);
        if (next_level.size() == 0) {
            auto seq_break = _check_from_chain_ends(xb-1, yb+1);
            if(seq_break != nullptr) { tree_.add_data(seq_break); }
        }
        
        if(current_pair.node->type() == SS_NodeData::SS_Type::SS_BULGE) {
            std::vector<SS_NodeDataOP> part_of_nway;
            std::vector<SS_NodeDataOP> not_part_of_nway;
            for(auto const & n : next_level) {
                if(n->type() == SS_NodeData::SS_Type::SS_BULGE ||
                   n->type() == SS_NodeData::SS_Type::SS_HAIRPIN) { part_of_nway.push_back(n); }
                else { not_part_of_nway.push_back(n); }
            }
            next_level = not_part_of_nway;
            //std::cout << xb << " " << yb << " " << part_of_nway.size() << " " << not_part_of_nway.size() << std::endl;
            
            if(part_of_nway.size() > 0) {
                auto ss_chains = current_pair.node->ss_chains();
                auto bounds = current_pair.node->bounds();
                
                for(auto const & n : next_level) {
                    for(int i = 0; i < n->ss_chains().size(); i++) {
                        if(n->ss_chains()[i]->sequence().size() > 0) {
                            ss_chains.push_back(n->ss_chains()[i]);
                        }
                    }
                    
                    if(n->type() == SS_NodeData::SS_Type::SS_BULGE) {
                        int nxb = _map_back_to_index(n->bound_side(0, SS_NodeData::Bound_Side::RIGHT))+1;
                        int nyb = _map_back_to_index(n->bound_side(1, SS_NodeData::Bound_Side::LEFT)) -1;
                        std::vector<SS_NodeDataOP> next_level_2 = _build_tree_level(nxb, nyb);
                        std::cout << next_level_2.size() << std::endl;
                        for(auto const & n2 : next_level_2) { next_level.push_back(n2); }
                    }
                }
                
                auto new_node = std::make_shared<SS_NodeData>(SS_NodeData::SS_Type::SS_NWAY,
                                                              ss_chains, bounds);
                index = current_pair.index;
                current_pair = NodeandIndex(new_node, index);
                
            }
        
        }
        
        index = tree_.add_data(current_pair.node ,current_pair.index);
        for (auto const & n : next_level) {
            open_nodes.push(NodeandIndex(n, index));
        }
        
    }
}
    
    
/*
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

*/

std::vector<SS_NodeDataOP>
SS_Tree::_build_tree_level(
    int xb,
    int yb) {
    
    std::vector<SS_NodeDataOP> next_level;
    while(xb <= yb) {
        if(xb == yb && residues_[xb]->dot_bracket() != ".") { break; }
        auto child = _assign_new_node(xb, yb);
        xb = _map_back_to_index(child->bound_side(1, SS_NodeData::Bound_Side::RIGHT))+1;
        next_level.push_back(child);
    }
    
    return next_level;
    
}
    
SS_NodeDataOP
SS_Tree::_assign_new_node(
    int xb,
    int yb) {
    
    std::vector<ResidueOPs> ss_chain_res(2);
    ResidueOPs bounds(2);
    bounds[0] = residues_.at(xb-1);
    bounds[1] = residues_.at(yb+1);
    int hairpin = 0;
    
    if(residues_.at(xb)->dot_bracket() == ".") {
        int end_x = _get_dot_bounds(xb, 0);
        for(int i = xb; i <= end_x; i++) {
            ss_chain_res[0].push_back(residues_[i]);
        }
        if(end_x == yb) { hairpin = 1; }
    }

    if(residues_.at(yb)->dot_bracket()== "." && hairpin == 0) {
        int end_y = _get_dot_bounds(yb, 1);
        for(int i = end_y; i <= yb; i++) {
            ss_chain_res[1].push_back(residues_[i]);
        }
    }
    
    SS_NodeData::SS_Type type;
    if(ss_chain_res[0].size() == 0 && ss_chain_res[1].size() == 0) {
        int pair = _get_brack_pair(xb);
        ss_chain_res[0].push_back(residues_.at(xb));
        ss_chain_res[1].push_back(residues_.at(pair));
        type = SS_NodeData::SS_Type::SS_BP;
        if(residues_.at(xb)->dot_bracket() == "[") {
            type = SS_NodeData::SS_Type::SS_PSEUDO_BP;
        }
    }
    
    else if(hairpin) {
        type = SS_NodeData::SS_Type::SS_HAIRPIN;
    }

    else {
        type = SS_NodeData::SS_Type::SS_BULGE;
    }

    ChainOPs ss_chains(2);
    ss_chains[0] = std::make_shared<Chain>(ss_chain_res[0]);
    ss_chains[1] = std::make_shared<Chain>(ss_chain_res[1]);
    
    auto current = std::make_shared<SS_NodeData>(type, ss_chains, bounds);
    
    return current;

}

    
int
SS_Tree::_get_brack_pair(
    int pos) {
    
    int bracket_count = 0;
    int i = -1;
    for(auto const & r : residues_) {
        i++;
        if( i < pos) { continue; }
        if (r->dot_bracket() == "(" || r->dot_bracket() == "[") { bracket_count++; }
        if (r->dot_bracket() == ")" || r->dot_bracket() == "]") { bracket_count--; }
        if(bracket_count == 0) { return i; }
    }
    
    throw SecondaryStructureException("could not find bracket pair");

}

int
SS_Tree::_get_dot_bounds(
    int pos,
    int reverse) {
    
    
    int incr = 1;
    if(reverse) { incr = -1; }
    
    while(1) {
        if(pos == 0 || pos == residues_.size()) { break; }
        if(residues_[pos]->dot_bracket() == ".") { pos += incr; }
        else{
            pos -= incr;
            break;
        }
    }
    return pos;
}
    
int
SS_Tree::_map_back_to_index(
    ResidueOP const & r) {
    
    int i = 0;
    try {
        i = (int)(std::find(residues_.begin(), residues_.end(), r) - residues_.begin());
    }
    catch(...) {
        throw SecondaryStructureException("unexpected error, cannot find residue in internal array, in SS_Tree, this should not happen!");
    }
    
    return i;
    
}
    
SS_NodeDataOP
SS_Tree::_check_from_chain_ends(
    int xb,
    int yb) {
    
    if(! _is_res_end_of_chain(residues_[xb]->num()) &&
       ! _is_res_end_of_chain(residues_[yb]->num())) {
        return nullptr;
    }
    
    auto fake_res_1 = std::make_shared<Residue>("&", "&", residues_[xb]->num(), "N", Uuid());
    auto fake_res_2 = std::make_shared<Residue>("&", "&", residues_[yb]->num(), "N", Uuid());
    auto ss_chain_res = std::vector<ResidueOPs>{ {fake_res_1}, {fake_res_2} };
    auto ss_chains = ChainOPs { std::make_shared<Chain>(ss_chain_res[0]),
                                std::make_shared<Chain>(ss_chain_res[1]) };
    auto bounds = ResidueOPs { residues_[xb-1], residues_[yb+1] };
    auto current = std::make_shared<SS_NodeData>(SS_NodeData::SS_Type::SS_SEQ_BREAK,
                                                 ss_chains, bounds);
    return current;
}

int
SS_Tree::_is_res_end_of_chain(
    int num) {
    
    for (auto const & c : ss_.chains()) {
        if(c->last()->num() == num) {
            return 0;
        }
    }
    return 1;
}

    
} //sstruct

