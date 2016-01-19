
//  motif_merger.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/4/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_data_structures/motif_merger.h"
#include "motif/motif_to_secondary_structure.h"

void
MotifMerger::add_motif(MotifOP const & m) {
    auto new_chains = ChainOPs(m->chains().size());
    int i =0;
    for(auto const & c : m->chains()) {
        new_chains[i] = std::make_shared<Chain>(*c);
        i++;
    }
    
    for(auto const & c : new_chains) {
        auto data = ChainNodeData(c, m->id());
        graph_.add_data(data, -1, -1, -1, 2, 1);
    }
    
    for(auto const & bp : m->basepairs()) {
        all_bps_[bp->uuid()] = bp;
    }
    
    motifs_[m->id()] = m;
    rebuild_structure_ = 1;
    
}

void
MotifMerger::add_motif(
    MotifOP const & m ,
    BasepairOP const & m_end,
    MotifOP const & parent,
    BasepairOP const & parent_end) {
    
    add_motif(m);
    _link_motifs(m, m_end, parent, parent_end);
    
    
}

void
MotifMerger::_link_motifs(
    MotifOP const & m1,
    BasepairOP const & m1_end,
    MotifOP const & m2,
    BasepairOP const & m2_end) {
    
    auto m1_end_nodes = _get_end_nodes(graph_.nodes(), m1_end);
    auto m2_end_nodes = _get_end_nodes(graph_.nodes(), m2_end);
    
    //std::cout << m2_end_nodes[0]->index() << " " << m2_end_nodes[1]->index() << std::endl;
    
    if(m2->mtype() == MotifType::HELIX && m1->mtype() != MotifType::HELIX) {
        _link_chains(m1_end_nodes, m2_end_nodes);
        bp_overrides_[m2_end->uuid()] = m1_end->uuid();
    }
    
    else {
        _link_chains(m2_end_nodes, m1_end_nodes);
        bp_overrides_[m1_end->uuid()] = m2_end->uuid();
    }

}

ChainNodes
MotifMerger::_get_end_nodes(
    ChainNodes const & nodes,
    BasepairOP const & end) {
    
    //auto end_nodes = std::make_unique<ChainNodes>(2);
    auto end_nodes = ChainNodes(2);
    for(auto const & n : nodes) {
        for(auto const & r : end->residues()) {
            if(n->data().c->first()->uuid() == r->uuid() && end_nodes[0] == nullptr) {
                end_nodes[0] = n;
            }
            else if(n->data().c->first()->uuid() == r->uuid() ) {
                throw MotifMergerException("cannot build chain map two residues are assigned \
                                           to 5' chain");
            }
            
            if(n->data().c->last()->uuid() == r->uuid() && end_nodes[1] == nullptr) {
                end_nodes[1] = n;
            }
            
            else if(n->data().c->last()->uuid() == r->uuid()) {
                throw MotifMergerException("cannot build chain map two residues are assigned \
                                           to 3' chain");
            }
        }
    }
    
    if(end_nodes[0] == nullptr || end_nodes[1] == nullptr) {
        throw MotifMergerException("did not build map properly, both chains are not found");
    }
    
    return end_nodes;
    
}

void
MotifMerger::_link_chains(
    ChainNodes & dominant_nodes,
    ChainNodes & auxiliary_nodes) {
    
    
    if(dominant_nodes[0] == dominant_nodes[1]) {
        _connect_chains(dominant_nodes[0], auxiliary_nodes[0], 1, 0);
        _connect_chains(dominant_nodes[0], auxiliary_nodes[1], 0, 1);
    }
    
    else if(auxiliary_nodes[0] == auxiliary_nodes[1]) {
        _connect_chains(dominant_nodes[1], auxiliary_nodes[0], 1, 0);
        _connect_chains(dominant_nodes[0], auxiliary_nodes[0], 0, 1);
    }
    
    else {
        _connect_chains(dominant_nodes[1], auxiliary_nodes[0], 1, 0);
        _connect_chains(dominant_nodes[0], auxiliary_nodes[1], 0, 1);
    }
}

void
MotifMerger::_connect_chains(
    ChainNode & d_node,
    ChainNode & a_node,
    int d_i,
    int a_i) {
    
    if(a_i == 0) {
        a_node->data().prime5_override = 1;
        res_overrides_[a_node->data().c->first()->uuid()] = d_node->data().c->last()->uuid();
    }
    else {
        a_node->data().prime3_override = 1;
        res_overrides_[a_node->data().c->last()->uuid()] = d_node->data().c->first()->uuid();
    }
    graph_.connect(d_node->index(), a_node->index(), d_i, a_i);
}


void
MotifMerger::_build_structure() {
    auto starts = ChainNodes();
    for(auto const & n : graph_.nodes()) {
        if(n->available_pos(0)) { starts.push_back(n); }
    }
    
    if(starts.size() == 0 and graph_.nodes().size() > 0) {
        throw MotifMergerException("no place to start in chain graph to build structure");
    }
    
    auto chains = ChainOPs();
    auto cur = ChainNode();
    for(auto const & n : starts) {
        auto res = ResidueOPs();
        cur = n;
        while(cur != nullptr) {
            auto new_res = cur->data().included_res();
            std::copy(new_res.begin(),
                      new_res.end(),
                      std::inserter(res, res.end()));
            if(cur->available_pos(1)) {
                cur = nullptr;
            }
            else {
                cur = cur->connections()[1]->partner(cur->index());
            }
        }
        chains.push_back(std::make_shared<Chain>(res));
    }
    
    auto struc = std::make_shared<Structure>(chains);
    auto all_bps_vec = BasepairOPs(all_bps_.size());
    int i = 0;
    for(auto const & kv : all_bps_) { all_bps_vec[i] = kv.second; i++; }
    auto bps = subselect_basepairs_with_res(struc->residues(), all_bps_vec);
    auto ends = end_from_basepairs(struc, *bps);
    rna_structure_ = std::make_shared<RNAStructure>(struc, *bps, *ends);
}

RNAStructureOP const &
MotifMerger::get_structure() {
    if(rebuild_structure_) {
        _build_structure();
        rebuild_structure_ = 0;
    }
    
    return rna_structure_;
}

void
MotifMerger::remove_motif(MotifOP const & m) {
    for(auto const & end : m->ends()) {
        if(bp_overrides_.find(end->uuid()) != bp_overrides_.end() ) {
            bp_overrides_.erase(end->uuid());
        }
        auto remove = std::vector<Uuid>();
        for(auto const & kv : bp_overrides_) {
            if(kv.second == end->uuid()) { remove.push_back(kv.first); }
        }
        for(auto const & u : remove) {
            bp_overrides_.erase(u);
        }
        
    }
    
    for(auto const & bp : m->basepairs()) {
        all_bps_.erase(bp->uuid());
    }
    
    auto remove = ChainNodes();
    for(auto const & n : graph_.nodes()) {
        if(n->data().m_id == m->id()) {
            remove.push_back(n);
        }
    }
    
    
    for(auto const & r : remove) {
        for(auto const & c : r->connections()) {
            if(c == nullptr) { continue; }
            auto p = c->partner(r->index());
            auto p_i = c->end_index(p->index());
            
            if(p_i == 0 && p->data().prime5_override == 1) {
                res_overrides_.erase(p->data().c->first()->uuid());
                p->data().prime5_override = 0;
            }
            else if(p_i == 1 && p->data().prime3_override == 1) {
                res_overrides_.erase(p->data().c->last()->uuid());
                p->data().prime3_override = 0;
            }
        }
        
        for(auto const & res: r->data().c->residues()) {
            if(res_overrides_.find(res->uuid()) != res_overrides_.end()) {
                res_overrides_.erase(res->uuid());
            }
            
        }
        
        graph_.remove_node(r->index());
    }
    
    auto it = motifs_.find(m->id());
    motifs_.erase(it);
    //motifs_.erase(m->id());
    rebuild_structure_ = 1;
    
}

void
MotifMerger::update_motif(MotifOP const & m) {
    
    for(auto const & n : graph_.nodes()) {
        if(n->data().m_id != m->id()) { continue; }
        auto res = n->data().c->residues();
        auto new_res = ResidueOPs();
        for(auto const & r : res) {
            auto new_r = m->get_residue(r->uuid());
            if(new_r == nullptr) {
                throw MotifMergerException("could not find res by uuid during update_motif");
            }
            new_res.push_back(new_r);
        }
        n->data().c = std::make_shared<Chain>(new_res);
    }
    
    for(auto const & bp : m->basepairs()) {
        all_bps_[bp->uuid()] = bp;
    }
    
    motifs_[m->id()] = m;
}


sstruct::PoseOP
MotifMerger::secondary_structure() {
    MotiftoSecondaryStructure parser;
    auto ss = parser.to_secondary_structure(get_structure());
    auto ss_motifs = sstruct::MotifOPs();
   
    auto r_cur = ResidueOP(nullptr);
    auto current_bp = BasepairOP(nullptr);
    auto ss_bp = sstruct::BasepairOPs();
    for(auto const & kv : motifs_) {
        auto m = kv.second;
        auto ss_chains = sstruct::ChainOPs();
        for(auto const & c : kv.second->chains()) {
            auto ss_res = sstruct::ResidueOPs();
            for(auto const & r : c->residues()) {
                r_cur = r;
                if(res_overrides_.find(r->uuid()) != res_overrides_.end()) {
                    r_cur = get_residue(res_overrides_[r->uuid()]);
                }
                auto ss_r = ss->get_residue(r_cur->uuid());
                if(ss_r == nullptr) {
                    throw MotifMergerException("could not find residue during ss build");
                }
                ss_res.push_back(ss_r);
            }
            
            ss_chains.push_back(std::make_shared<sstruct::Chain>(ss_res));
        }
        
        auto ss_struc = std::make_shared<sstruct::Structure>(ss_chains);
        auto ss_bps = sstruct::BasepairOPs();
        for(auto const & bp : kv.second->basepairs()) {
            if(bp->bp_type() != "cW-W")     { continue; }
            if(! wc_bp(bp) && ! gu_bp(bp) ) { continue; }
            current_bp = bp;
            if(bp_overrides_.find(current_bp->uuid()) != bp_overrides_.end()) {
                current_bp = get_basepair(bp_overrides_[current_bp->uuid()]);
            }
            ss_bp = ss->get_basepair(current_bp->uuid());
            if(ss_bp.size() == 0) {
                throw MotifMergerException("could not find basepair during ss build");
            }
            ss_bps.push_back(ss_bp[0]);
        }
        
        auto ss_ends = sstruct::BasepairOPs();
        for(auto const & end : kv.second->ends()) {
            current_bp = end;
            if(bp_overrides_.find(current_bp->uuid()) != bp_overrides_.end()) {
                current_bp = get_basepair(bp_overrides_[current_bp->uuid()]);
            }
            ss_bp = ss->get_basepair(current_bp->uuid());
            if(ss_bp.size() == 0) {
                throw MotifMergerException("could not find basepair during ss build");
            }
            ss_ends.push_back(ss_bp[0]);
        }
        
        auto ss_motif = std::make_shared<sstruct::Motif>(ss_struc, ss_bps, ss_ends, m->end_ids(),
                                                         m->name(), m->path(), m->score());
        ss_motif->id(m->id());
        ss_motif->mtype(m->mtype());
        ss_motifs.push_back(ss_motif);
    }
    
    return std::make_shared<sstruct::Pose>(ss, ss_motifs);
}
































