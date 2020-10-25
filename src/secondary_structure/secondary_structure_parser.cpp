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

        if(n->data().residues[0]->dot_bracket() != "(") {
            continue;
        }

        auto nodes = std::vector<SSNodeOP>();
        auto type = _walk_nodes_new(n, nodes);
        for(auto const & n : nodes) { seen_[n] = 1; }

        if(type == util::MotifType::UNKNOWN) { continue; }
        auto s = _generate_structure_from_nodes(nodes);
        auto m = _build_motif(s);
        if(m == nullptr) { continue; }
        motifs.push_back(m);
        //std::cout << m->sequence() << " " << m->dot_bracket() << std::endl;

        /*auto m = _generate_motif(n);
        seen_[n] = 1;
        if(m == nullptr) { continue; }
        motifs.push_back(m);*/
    }

    for(auto const & n : *g) {
        if(seen_.find(n) == seen_.end()) {
            auto s = _generate_structure_from_nodes(std::vector<SSNodeOP>{n});
            auto m = _build_motif(s);
            motifs.push_back(m);
        }
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


StructureOP
Parser::_generate_structure_from_nodes(
        std::vector<SSNodeOP> const & nodes) {
    auto res = ResidueOPs();
    for(auto const & n : nodes) {
        for(auto const & r : n->data().residues) {
            res.push_back(r);
        }
    }
    auto chains = ChainOPs();
    auto current_chain_res = ResidueOPs();
    current_chain_res.push_back(res[0]);
    auto diff = 0;
    for(int i = 1; i < res.size(); i++) {
        diff = res[i]->num() - res[i-1]->num();
        if(diff != 1) {
            chains.push_back(std::make_shared<Chain>(current_chain_res));
            current_chain_res = ResidueOPs();
            current_chain_res.push_back(res[i]);
        }
        else {
            current_chain_res.push_back(res[i]);
        }
    }
    chains.push_back(std::make_shared<Chain>(current_chain_res));

    return std::make_shared<Structure>(chains);

}

util::MotifType
Parser::_walk_nodes_new(
        SSNodeOP const & n,
        std::vector<SSNodeOP> & nodes) {

    int i = 0;
    auto current = n;
    auto next = SSNodeOP(nullptr);

    auto current_bp_res = ResidueOPs(2);
    auto next_bp_res = ResidueOPs(2);
    while(current != nullptr) {
        nodes.push_back(current);
        if(current->connections()[1] != nullptr) {
            next = current->connections()[1]->partner(current->index());
        }
        else {
            break;
        }

        /*std::cout << "CURRENT: " << current->data().residues[0]->num() << " NEXT: " << next->data().residues[0]->num() <<  " LIST: ";
        for(auto const & n: nodes) {
            std::cout << n->data().residues[0]->num() << " ";
        }
        std::cout << std::endl;
        */

        if(current->data().type == NodeType::PAIRED && next->data().type == NodeType::PAIRED) {
            _get_basepair_res(current, current_bp_res);
            _get_basepair_res(next, next_bp_res);

            // check if base pair step
            if (current_bp_res[0]->num() == next_bp_res[0]->num() - 1 &&
                current_bp_res[1]->num() == next_bp_res[1]->num() + 1) {
                nodes.push_back(next);
                nodes.push_back(next->connections()[2]->partner(next->index()));
                nodes.push_back(current->connections()[2]->partner(current->index()));
                return util::MotifType::HELIX;
            }

            nodes.push_back(next);

            if(_is_a_bp(nodes[0], next)) {
                return util::MotifType::NWAY;
            }

            current = next;
            next = next->connections()[2]->partner(next->index());

        }

        else if(current->data().type == NodeType::PAIRED && next->data().type == NodeType::UNPAIRED) {
            if(std::find(nodes.begin(), nodes.end(), next) != nodes.end()) {
                return util::MotifType::UNKNOWN;
            }
        }

        else if(current->data().type == NodeType::UNPAIRED && next->data().type == NodeType::PAIRED) {
            // come full circle
            if(_is_a_bp(nodes[0], next)) {
                nodes.push_back(next);
                if(nodes.size() == 3)      { return util::MotifType::HAIRPIN; }
                else if(nodes.size() == 5) { return util::MotifType::TWOWAY; }
                else                       { return util::MotifType::NWAY; }
            }

            if(std::find(nodes.begin(), nodes.end(), next) != nodes.end()) {
                return util::MotifType::UNKNOWN;
            }

            nodes.push_back(next);
            current = next;
            next = next->connections()[2]->partner(next->index());

        }


        current = next;

    }
    return util::MotifType::UNKNOWN;
}

void
Parser::_get_basepair_res(
        SSNodeOP n,
        ResidueOPs & bp_res) {
    if(n->data().type == NodeType::UNPAIRED) {
        throw secondary_structure::Exception("cannot get baepair res it is not a basepair");
    }
    auto partner =n->connections()[2]->partner(n->index());
    bp_res[0] = n->data().residues[0];
    bp_res[1] = partner->data().residues[0];

}


bool
Parser::_is_a_bp(
        SSNodeOP n1,
        SSNodeOP n2) {
    if(n1->connections()[2] == nullptr) { return false; }
    if(n1->connections()[2]->partner(n1->index()) == n2) { return true;}
    return false;
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
    if(last_node->connections()[2] == nullptr) {
        return nullptr;
    }

    return last_node->connections()[2]->partner(last_node->index());
}

    
MotifOP
Parser::_generate_motif(
    SSNodeOP const & n) {
    
    auto next_n = _walk_nodes(n);
    auto chains = ChainOPs{chain_};

    int i =0 ;
    if(next_n != nullptr) {
        while (next_n != n) {
            next_n = _walk_nodes(next_n);
            if (next_n == nullptr) { return nullptr; }
            chains.push_back(chain_);
            std::cout << chains.size() << " " << next_n->data().residues[0]->num() << std::endl;
            i += 1;
            if(i > 10) {
                exit(0);
            }
        }
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
    if(ends.size() > 0) {
        auto end_ids = Strings();
        for (auto const & end : m->ends()) {
            auto end_id = assign_end_id(m, end);
            end_ids.push_back(end_id);
        }
        m->end_ids(end_ids);
    }
    
    if(m->residues().size() == 4) {
        m->mtype(util::MotifType::HELIX);
    }
    else if(m->chains().size() == 2) {
        m->mtype(util::MotifType::TWOWAY);
    }
    else if(m->chains().size() == 1 && m->basepairs().size() == 0) {
        m->mtype(util::MotifType::SSTRAND);
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





































