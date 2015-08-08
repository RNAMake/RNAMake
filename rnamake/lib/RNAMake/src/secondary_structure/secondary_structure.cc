//
//  secondary_structure.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <memory.h>
#include <map>
#include <queue>
#include <sstream>

#include "secondary_structure/chain.h"
#include "secondary_structure/secondary_structure.h"

namespace sstruct {


SecondaryStructure::SecondaryStructure(
    String const & sequence,
    String const & dot_bracket):
    Motif()
    {
    
    if(sequence.length() != dot_bracket.length()) {
        throw SecondaryStructureException("cannot construct new SecondaryStructure object: new sequence and dot bracket are not the same length");
    }

    if(sequence.length() == 0) {
        throw SecondaryStructureException("cannot construct new SecondaryStructure object: sequence is of lenght zero!");
    }
    
    if(dot_bracket[0] != '(' && dot_bracket[0] != '.' &&
       dot_bracket[0] != '&' && dot_bracket[0] != '+') {
        throw SecondaryStructureException("cannot construct new SecondaryStructure object: dot bracket notation for secondary structure is not valid. perhaps you flipped seq and ss?");
    }
    
    ResidueOPs res;
    int count = 0;
    String chain_ids = "ABCDEFGHIKLMNO";
    int i = -1, ci = 0;
    for(auto & s : sequence) {
        i++;
        if(s != '&' && s != '+') {
            String name = "", db = "", chain_id = "";
            name += s; db += dot_bracket[i]; chain_id += chain_ids[ci];
            auto r = std::make_shared<Residue>(name, db, count, chain_id, Uuid());
            res.push_back(r);
            count++;
        }
        else {
            ci += 1;
            this->chains_.push_back(std::make_shared<Chain>(res));
            res = ResidueOPs();
        }
    }
    
    if(res.size() > 0) {
        this->chains_.push_back(std::make_shared<Chain>(res));
    }
    
}

SecondaryStructure::SecondaryStructure(
    ChainOPs const & chains) {
 
    this->chains_ = chains;

}

String
SecondaryStructure::to_str() const {
    std::stringstream ss;
    for(auto const & c : chains_) {
        ss << c->to_str() + ":";
    }
    ss << "@";
    for(auto const & bp : basepairs_) {
        ss << bp->res1()->num() << "|" << bp->res1()->chain_id() << "|" << bp->res1()->i_code() << "_";
        ss << bp->res2()->num() << "|" << bp->res2()->chain_id() << "|" << bp->res2()->i_code() << ",";
    }
    ss << "@";
    for(auto const & bp : ends_) {
        ss << bp->res1()->num() << "|" << bp->res1()->chain_id() << "|" << bp->res1()->i_code() << "_";
        ss << bp->res2()->num() << "|" << bp->res2()->chain_id() << "|" << bp->res2()->i_code() << ",";
    }
    ss << "@";
    for(auto const & kv : motifs_) {
        if(kv.first == "ALL") { continue; }
        for(auto const & m : kv.second) {
            ss << m->to_str() + "$";
        }
    }
    
    return ss.str();
    
}
    
String
assign_end_id(
    SecondaryStructureOP const & ss,
    BasepairOP const & end) {
    
    int found = 0;
    for(auto const & e : ss->ends()) {
        if(end == e) { found = 1; break; }
    }
    
    if(!found) {
        throw SecondaryStructureException("supplied an end that is not in current ss element");
    }
    
    ChainOPs all_chains;
    for(auto const & c : ss->chains()) { all_chains.push_back(c); }
    std::queue<ChainOP> open_chains;
    
    for(auto const & c : all_chains) {
        if(c->first() == end->res1() || c->first() == end->res2()) {
            open_chains.push(c);
            break;
        }
    }
    
    all_chains.erase(std::remove(all_chains.begin(), all_chains.end(), open_chains.front()),
                     all_chains.end());
    
    std::map<ResidueOP, int> seen_res;
    std::map<BasepairOP, int> seen_bp;
    BasepairOP saved_bp, bp;
    String structure, sequence, dot_bracket;
    using IdChain = std::pair<String, String>;
    std::vector<IdChain> id_chains;
    int score = 0, pos = 0, i =0;
    ChainOPs best_chains;
    ChainOP best_chain;
    while( ! open_chains.empty() ) {
        ChainOP c = open_chains.front();
        open_chains.pop();
        
        for(auto const & r : c->residues()) {
            dot_bracket = ".";
            bp = ss->get_bp(r);
            saved_bp = nullptr;
            if(bp != nullptr) {
                saved_bp = bp;
                auto partner_r = bp->partner(r);
                if     (seen_bp.find(bp) == seen_bp.end() &&
                        seen_res.find(r) == seen_res.end() &&
                        seen_res.find(partner_r) == seen_res.end()) {
                    seen_res[r] = 1;
                    dot_bracket = "(";
                }
                else if(seen_res.find(partner_r) != seen_res.end()) {
                    if(seen_res[r] > 1) {
                        dot_bracket = ".";
                    }
                    else {
                        dot_bracket = ")";
                        seen_res[r] = 1;
                        seen_res[partner_r] += 1;
                    }
                }
            }
            structure += dot_bracket;
            sequence  += r->name();
            
            if(saved_bp != nullptr) { seen_bp[saved_bp] = 1; }
        }
        
        id_chains.push_back(IdChain(sequence, structure));
        sequence = ""; structure = "";
        int best_score = -1;
        
        for(auto const & c : all_chains) {
            score = 0;
            for(auto const & r : c->residues()) {
                bp = ss->get_bp(r);
                if(bp != nullptr && seen_bp.find(bp) != seen_bp.end()) {
                    score += 1;
                }
            }
            if(score > best_score) { best_score = score; }
        }
        
        best_chains.resize(0);
        for(auto const & c : all_chains) {
            score = 0;
            for(auto const & r : c->residues()) {
                bp = ss->get_bp(r);
                if(bp != nullptr && seen_bp.find(bp) != seen_bp.end()) {
                    score += 1;
                }
            }
            if (score == best_score) {
                best_chains.push_back(c);
            }
        }
        
        best_score = 10000;
        best_chain = nullptr;
        for(auto const & c : best_chains) {
            pos = 1000;
            i = 0;
            for(auto const & r : c->residues()) {
                bp = ss->get_bp(r);
                if(bp != nullptr && seen_bp.find(bp) != seen_bp.end()) {
                    pos = i;
                    break;
                }
                i++;
            }
            if(pos < best_score) {
                best_score = pos;
                best_chain = c;
            }
        }
        
        if(best_chain == nullptr) { break; }
        all_chains.erase(std::remove(all_chains.begin(), all_chains.end(), best_chain),
                         all_chains.end());
        open_chains.push(best_chain);
    }
    
    String ss_id = "";
    i = 0;
    for(auto const & id_chain : id_chains) {
        ss_id += id_chain.first + "_";
        for(auto const & e : id_chain.second) {
            if     (e == '(') {
                ss_id += "L";
            }
            else if(e == ')') {
                ss_id += "R";
            }
            else if(e == '.') {
                ss_id += ".";
            }
            else {
                throw SecondaryStructureException("unexpected symbol in dot bracket notation: " + e);
            }
        }
        if(i != id_chains.size()-1) { ss_id += "_"; }
        i++;
    }
    
    return ss_id;
    
}

    
ResidueOP
get_res_from_res_str(
    SecondaryStructure const & ss,
    String const & res_str) {
    
    Strings res_spl = split_str_by_delimiter(res_str, "|");
    int num = std::stoi(res_spl[0]);
    if(res_str.size() == 3) {
        return ss.get_residue(num, res_spl[1], res_spl[2]);
    }
    else {
        return ss.get_residue(num, res_spl[1]);

    }
}
    

SecondaryStructure
str_to_secondary_structure(
    String const & s) {
    
  
    auto spl = split_str_by_delimiter(s, "@");
    auto c_spl = split_str_by_delimiter(spl[0], ":");
    ChainOPs chains;
    for(auto const & chain_str : c_spl) {
        auto c = std::make_shared<Chain>(str_to_chain(chain_str));
        chains.push_back(c);
    }
    auto ss = SecondaryStructure(chains);
    auto bp_spl = split_str_by_delimiter(spl[1], ",");
    BasepairOPs basepairs, ends;
    ResidueOP res_1, res_2;

    for(auto const & bp_str : bp_spl) {
        auto res_info = split_str_by_delimiter(bp_str, "_");
        res_1 = get_res_from_res_str(ss, res_info[0]);
        res_2 = get_res_from_res_str(ss, res_info[1]);
        auto bp = std::make_shared<Basepair>(res_1, res_2, Uuid());
        basepairs.push_back(bp);
    }
    ss.basepairs(basepairs);
    auto end_spl = split_str_by_delimiter(spl[2], ",");
    for(auto const & end_str : end_spl) {
        auto res_info = split_str_by_delimiter(end_str, "_");
        res_1 = get_res_from_res_str(ss, res_info[0]);
        res_2 = get_res_from_res_str(ss, res_info[1]);
        auto bp = ss.get_bp(res_1, res_2);
    }
    if(spl.size() == 3) {
        return ss;
    }
    
    auto motif_spl = split_str_by_delimiter(spl[3], "$");
    std::map<String, MotifOPs> motifs;
    motifs["ALL"] = MotifOPs();
    for(auto const & motif_str : motif_spl) {
        auto motif_info = split_str_by_delimiter(motif_str, ";");
        auto type = motif_info[0];
        auto chain_spl = split_str_by_delimiter(motif_info[1], ",");
        ChainOPs m_chains;
        for(auto const & res_list : chain_spl) {
            auto res = split_str_by_delimiter(res_list, "_");
            ResidueOPs residues;
            for(auto const & r_str : res) {
                auto r = get_res_from_res_str(ss, r_str);
                residues.push_back(r);
            }
            m_chains.push_back(std::make_shared<Chain>(residues));
        }
        BasepairOPs m_basepairs, m_ends;
        auto bp_spl = split_str_by_delimiter(motif_info[2], ",");
        for(auto const & bp_str : bp_spl) {
            auto res_info = split_str_by_delimiter(bp_str, "_");
            res_1 = get_res_from_res_str(ss, res_info[0]);
            res_2 = get_res_from_res_str(ss, res_info[1]);
            auto bp = ss.get_bp(res_1, res_2);
            if(bp == nullptr) {
                throw std::runtime_error("cannot find basepair in str_to_ss");
            }
            m_basepairs.push_back(bp);
        }
        auto end_spl = split_str_by_delimiter(motif_info[3], ",");
        for(auto const & end_str : end_spl) {
            auto res_info = split_str_by_delimiter(end_str, "_");
            res_1 = get_res_from_res_str(ss, res_info[0]);
            res_2 = get_res_from_res_str(ss, res_info[1]);
            auto bp = ss.get_bp(res_1, res_2);
            if(bp == nullptr) {
                throw std::runtime_error("cannot find basepair in str_to_ss");
            }
            m_ends.push_back(bp);
        }
        auto m = std::make_shared<Motif>(type, m_ends, m_chains);
        m->basepairs(m_basepairs);
        motifs["ALL"].push_back(m);
        if(motifs.find(m->type()) == motifs.end()) {
            motifs[m->type()] = MotifOPs();
        }
        motifs[m->type()].push_back(m);
    }
    
    ss.set_motifs(motifs);
    
    return ss;
    
}
    
}






















