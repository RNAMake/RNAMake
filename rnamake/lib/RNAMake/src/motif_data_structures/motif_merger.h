//
//  motif_merger.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/4/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_merger__
#define __RNAMake__motif_merger__

#include <stdio.h>
#include <map>

//RNAMake Headers
#include "data_structure/graph/graph.h"
#include "util/uuid.h"
#include "secondary_structure/pose.h"
#include "structure/residue.h"
#include "structure/chain.fwd.h"
#include "structure/chain.h"
#include "structure/rna_structure.h"
#include "motif/motif.h"


class MotifMergerException : public std::runtime_error {
public:
    MotifMergerException(
        String const & message) :
    std::runtime_error(message)
    {}
};

struct ChainNodeData {
    inline
    ChainNodeData() {}
    
    inline
    ChainNodeData(
        ChainOP const & nc,
        Uuid const & nm_id,
        int nprime5_override = 0,
        int nprime3_override = 0):
    c(nc),
    m_id(nm_id),
    prime5_override(nprime5_override),
    prime3_override(nprime3_override)
    {}

    inline
    ChainNodeData(
        ChainNodeData const & cnd):
    c(cnd.c),
    m_id(cnd.m_id),
    prime5_override(cnd.prime5_override),
    prime3_override(cnd.prime3_override)
    {}

    ResidueOPs
    included_res() {
        auto res = c->residues();
        if(prime5_override == 1) { res.erase( res.begin() ); }
        if(prime3_override == 1) { res.pop_back(); }
        return res;
    }
    
    ChainOP c;
    Uuid m_id;
    int prime5_override, prime3_override;
};

typedef GraphNodeOP<ChainNodeData>   ChainNode;
typedef std::vector<ChainNode>       ChainNodes;
typedef std::unique_ptr<ChainNodes>  ChainNodesUP;

class MotifMerger {
public:
    MotifMerger():
    all_bps_(std::map<Uuid, BasepairOP, UuidCompare>()),
    motifs_(std::map<Uuid, MotifOP, UuidCompare>()),
    res_overrides_(std::map<Uuid, Uuid, UuidCompare>()),
    bp_overrides_(std::map<Uuid, Uuid, UuidCompare>()),
    graph_(GraphStatic<ChainNodeData>()),
    rebuild_structure_(1),
    rna_structure_(std::make_shared<RNAStructure>()){
        
    }
    
    MotifMerger(
        MotifMerger const & mm,
        MotifOPs const & motifs):
    all_bps_(std::map<Uuid, BasepairOP, UuidCompare>()),
    motifs_(std::map<Uuid, MotifOP, UuidCompare>()),
    res_overrides_(mm.res_overrides_),
    bp_overrides_(mm.bp_overrides_),
    rebuild_structure_(1) {
        
        graph_ = GraphStatic<ChainNodeData>(mm.graph_);
        res_overrides_ = mm.res_overrides_;
        for(auto const & n : mm.graph_.nodes()) {
            graph_.get_node(n->index())->data() = ChainNodeData(n->data());
  
        }
        
        for(auto const & m : motifs) { update_motif(m); }
    }
    
    ~MotifMerger() { }
    
public:
    void
    add_motif(MotifOP const &);
    
    void
    add_motif(
        MotifOP const &,
        BasepairOP const &,
        MotifOP const &,
        BasepairOP const &);
    
    RNAStructureOP const &
    get_structure();
    
    
    void
    remove_motif(MotifOP const &);
    
    void
    update_motif(MotifOP const &);
    
    void
    connect_motifs(
        MotifOP const & m1,
        MotifOP const & m2,
        BasepairOP const & m1_end,
        BasepairOP const & m2_end) {
        _link_motifs(m1, m1_end, m2, m2_end);
    }
    
    sstruct::PoseOP
    secondary_structure();
    
public:
    
    inline
    ResidueOP
    get_residue(Uuid const & uuid) {
        auto r = ResidueOP(nullptr);
        for(auto const & kv : motifs_) {
            r = kv.second->get_residue(uuid);
            if(r != nullptr) { return r; }
        }
        return nullptr;
    }
    
    inline
    BasepairOP
    get_basepair(Uuid const & uuid) {
        if(all_bps_.find(uuid) != all_bps_.end()) {
            return all_bps_[uuid];
        }
        else{
            return nullptr;
        }
    }
    
    inline
    void
    to_pdb(
        String const & fname,
        int renumber = -1) {
        return get_structure()->to_pdb(fname, renumber);
    }
    
private:
    void
    _link_motifs(
        MotifOP const &,
        BasepairOP const &,
        MotifOP const &,
        BasepairOP const &);
    
    ChainNodes
    _get_end_nodes(
        ChainNodes const &,
        BasepairOP const &);
    
    void
    _link_chains(
        ChainNodes &,
        ChainNodes &);
    
    void
    _connect_chains(
        ChainNode &,
        ChainNode &,
        int d_i,
        int a_i);
    
    void
    _build_structure();
    
private:
    std::map<Uuid, BasepairOP, UuidCompare> all_bps_;
    std::map<Uuid, MotifOP, UuidCompare> motifs_;
    std::map<Uuid, Uuid, UuidCompare> res_overrides_, bp_overrides_;
    GraphStatic<ChainNodeData> graph_;
    int rebuild_structure_;
    RNAStructureOP rna_structure_;
    
};

typedef std::shared_ptr<MotifMerger> MotifMergerOP;

#endif /* defined(__RNAMake__motif_merger__) */



































