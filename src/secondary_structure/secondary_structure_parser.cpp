//
//  secondary_structure_parser.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure/secondary_structure_parser.h"
#include "secondary_structure/util.h"

namespace secondary_structure {

SecondaryStructureChainGraphOP
Parser::parse(
    String const & sequence,
    String const & dot_bracket) {

    structure_ = std::make_shared<Structure>(sequence, dot_bracket);
    residues_  = structure_->residues();
    
    auto g = std::make_shared<SecondaryStructureChainGraph>();
    
    auto res = ResidueOPs();
    pairs_ = BasepairOPs();
    for(auto const & r : residues_) {
        int is_start_res = _start_of_chain(r);
        if(r->dot_bracket() == ".") { res.push_back(r); }
        
        else if(r->dot_bracket() == "(") {
            if(res.size() > 0) {
                _add_unpaired_residues_to_graph(g, res, is_start_res);
                res = ResidueOPs();
            }
           
            _add_paired_res_to_graph(g, r, is_start_res);
        }
        
        else if(r->dot_bracket() == ")") {
            if(res.size() > 0) {
                _add_unpaired_residues_to_graph(g, res, is_start_res);
                res = ResidueOPs();
            }
            
            auto pair = _get_previous_pair(r);
            auto new_data = NodeData(ResidueOPs{r}, NodeType::PAIRED);
            int parent_index = g->get_node_by_res(_previous_res(r));
            int pair_res_pos = g->get_node_by_res(pair->res1());
            int pos = g->add_chain(new_data, parent_index, is_start_res);
            g->pair_res(pair_res_pos, pos);
        }
        
        else {
            throw Exception("unexpected symbol in dot bracket: " + r->dot_bracket());
        }
    }
    
    if(res.size() > 0) { _add_unpaired_residues_to_graph(g, res, 0); }
    
    return g;
}
    
void
Parser::_add_unpaired_residues_to_graph(
    SecondaryStructureChainGraphOP & g,
    ResidueOPs const & res,
    int is_start_res) {
    
    auto parent_index = g->get_node_by_res(_previous_res(res[0]));
    auto new_data = NodeData(res, NodeType::UNPAIRED);
    g->add_chain(new_data, parent_index, is_start_res);

}

    
void
Parser::_add_paired_res_to_graph(
    SecondaryStructureChainGraphOP & g,
    ResidueOP const & r,
    int is_start_res) {
    
    auto pair_res = _get_bracket_pair(r);
    int parent_index = g->get_node_by_res(_previous_res(r));
    auto pair = std::make_shared<Basepair>(r, pair_res, util::Uuid());
    pairs_.push_back(pair);
    auto new_data = NodeData(ResidueOPs{r}, NodeType::PAIRED);
    g->add_chain(new_data, parent_index, is_start_res);
    
}
    
    
BasepairOP
Parser::_get_previous_pair(
    ResidueOP const & r)  {
    
    BasepairOP pair = nullptr;
    for(auto const & p : pairs_) {
        if(p->res2()->uuid() == r->uuid()) {
            pair = p;
            break;
        }
    }

    if(pair == nullptr) {
        throw Exception(
            "cannot parse secondary structure: \n" + structure_->sequence() + "\n" +
            structure_->dot_bracket() + "\n position: " + std::to_string(r->num()) +
            " has no matching pair");
    }
    
    return pair;
    
}


MotifOPs
Parser::parse_to_motifs(
    String const & sequence,
    String const & dot_bracket) {
    
    auto g = parse(sequence, dot_bracket);
    return _parse_to_motifs(g);
}

MotifOPs
Parser::_parse_to_motifs(
        SecondaryStructureChainGraphOP g) {
    auto motifs = MotifOPs();
    seen_ = std::map<SSNodeOP, int>();

    for(auto const & n : *g) {
        if(n->connections()[1] == nullptr) {
            continue;
        }
        if(n->data().type == NodeType::UNPAIRED) {
            continue;
        }
        auto m = _generate_motif(n);
        seen_[n] = 1;
        if(m == nullptr) { continue; }
        motifs.push_back(m);
    }
    return motifs;

}


MotifOP
Parser::parse_to_motif(
    String const & sequence,
    String const & dot_bracket) {
    
    parse(sequence, dot_bracket);
    return _build_motif(structure_);
    
}
    
PoseOP
Parser::parse_to_pose(
    String const & sequence,
    String const & dot_bracket) {
    
    auto motifs = parse_to_motifs(sequence, dot_bracket);
    auto m = _build_motif(structure_);

    return std::make_shared<Pose>(m, motifs);
}

    
SSNodeOP
Parser::_walk_nodes(
    SSNodeOP const & n) {
    
    int bps_count = 0;
    auto res = ResidueOPs();
    auto current = n;
    auto last_node = n;
    if(n == nullptr) { return nullptr; }
    while(current != nullptr) {
        if(seen_.find(current) != seen_.end()) {
            return nullptr;
        }
        if(current->data().type == NodeType::PAIRED) { bps_count += 1; }
        for(auto & r : current->data().residues) {
            res.push_back(r);
        }
        /*std::copy(current->data().residues.begin(),
                  current->data().residues.end(),
                  std::inserter(res, res.end()));*/
        last_node = current;
        if(current->connections()[1] != nullptr) {
            current = current->connections()[1]->partner(current->index());
        }
        else {
            break;
        }
        
        if(bps_count == 2) { break; }

    }
    chain_ = std::make_shared<Chain>(res);
    return last_node->connections()[2]->partner(last_node->index());
}

    
MotifOP
Parser::_generate_motif(
    SSNodeOP const & n) {
    
    auto next_n = _walk_nodes(n);
    auto chains = ChainOPs{chain_};
    
    while(next_n != n) {
        next_n = _walk_nodes(next_n);
        if(next_n == nullptr) { return nullptr; }
        chains.push_back(chain_);
    }
    
    auto struc = std::make_shared<Structure>(chains);
    if(struc->residues().size() < 3) {
        return nullptr;
    }    
 
    return _build_motif(struc);
}
    
MotifOP
Parser::_build_motif(
    StructureOP const & struc) {
    
    auto res = std::map<ResidueOP, int>();
    for(auto const & r : struc->residues()) {
        res[r] = 1;
    }
    
    auto bps = BasepairOPs();
    for(auto const & bp : pairs_) {
        if(res.find(bp->res1()) != res.end() &&
           res.find(bp->res2()) != res.end()) {
            bps.push_back(bp);
        }
    }
    
    auto chain_ends = std::map<ResidueOP, int>();
    for(auto const & c : struc->chains()) {
        chain_ends[c->first()] = 1;
        chain_ends[c->last()] = 1;
    }
    
    auto ends = BasepairOPs();
    for(auto const & bp : bps) {
        if(chain_ends.find(bp->res1()) != chain_ends.end() &&
           chain_ends.find(bp->res2()) != chain_ends.end()) {
            ends.push_back(bp);
        }
    }
    
    auto m = std::make_shared<Motif>(struc, bps, ends);
    auto end_ids = Strings();
    for(auto const & end : m->ends()) {
        auto end_id = assign_end_id(m, end);
        end_ids.push_back(end_id);
    }
    m->end_ids(end_ids);
    
    if(m->residues().size() == 4) {
        m->mtype(util::MotifType::HELIX);
    }
    else if(m->chains().size() == 2) {
        m->mtype(util::MotifType::TWOWAY);
    }
    else if(m->chains().size() == 1) {
        m->mtype(util::MotifType::HAIRPIN);
    }
    else {
        m->mtype(util::MotifType::NWAY);
    }
    
    
    return m;
    
}
    
}





































