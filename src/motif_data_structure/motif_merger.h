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

namespace motif_data_structure {

class MotifMergerException : public std::runtime_error {
public:
    MotifMergerException(
            String const & message) :
            std::runtime_error(message) {}
};


struct ChainNodeData {
    inline
    ChainNodeData() {}

    inline
    ChainNodeData(
            structure::ChainOP const & nc,
            util::Uuid const & nm_id,
            int nprime5_override = 0,
            int nprime3_override = 0) :
            c(nc),
            m_id(nm_id),
            prime5_override(nprime5_override),
            prime3_override(nprime3_override) {}

    inline
    ChainNodeData(
            ChainNodeData const & cnd) :
            c(cnd.c),
            m_id(cnd.m_id),
            prime5_override(cnd.prime5_override),
            prime3_override(cnd.prime3_override) {}

    structure::ResidueOPs
    included_res() {
        auto res = c->residues();
        if (prime5_override == 1) { res.erase(res.begin()); }
        if (prime3_override == 1) { res.pop_back(); }
        return res;
    }

    structure::ChainOP c;
    util::Uuid m_id;
    int prime5_override, prime3_override;
};

typedef data_structure::graph::GraphNodeOP<ChainNodeData> ChainNode;
typedef std::vector<ChainNode> ChainNodes;
typedef std::unique_ptr<ChainNodes> ChainNodesUP;

class MotifMerger {
public:
    MotifMerger() :
            all_bps_(std::map<util::Uuid, structure::BasepairOP, util::UuidCompare>()),
            bp_overrides_(std::map<util::Uuid, util::Uuid, util::UuidCompare>()),
            graph_(data_structure::graph::GraphStatic<ChainNodeData>()),
            motifs_(std::map<util::Uuid, motif::MotifOP, util::UuidCompare>()),
            rebuild_structure_(1),
            res_overrides_(std::map<util::Uuid, util::Uuid, util::UuidCompare>()),
            rna_structure_(std::make_shared<structure::RNAStructure>())
    {

    }

    MotifMerger(
            MotifMerger const & mm,
            motif::MotifOPs const & motifs) :
            all_bps_(std::map<util::Uuid, structure::BasepairOP, util::UuidCompare>()),
            motifs_(std::map<util::Uuid, motif::MotifOP, util::UuidCompare>()),
            res_overrides_(mm.res_overrides_),
            bp_overrides_(mm.bp_overrides_),
            rna_structure_(std::make_shared<structure::RNAStructure>()),
            graph_(data_structure::graph::GraphStatic<ChainNodeData>(mm.graph_)),
            rebuild_structure_(1) {

        for (auto const & n : mm.graph_.nodes()) {
            graph_.get_node(n->index())->data() = ChainNodeData(n->data());

        }

        for (auto const & m : motifs) { update_motif(m); }
    }

    ~MotifMerger() {}

private:
    enum MotifMergerType {
        // sequence matters for these motifs
                SPECIFIC_SEQUENCE = 0,
        // sequence does not matter such as idealized helices
                NON_SPECIFIC_SEQUENCE = 1
    };

public:
    void
    add_motif(motif::MotifOP const &);

    void
    add_motif(
            motif::MotifOP const &,
            structure::BasepairOP const &,
            motif::MotifOP const &,
            structure::BasepairOP const &);

    structure::RNAStructureOP const &
    get_structure();


    void
    remove_motif(motif::MotifOP const &);

    void
    update_motif(motif::MotifOP const &);

    void
    connect_motifs(
            motif::MotifOP const & m1,
            motif::MotifOP const & m2,
            structure::BasepairOP const & m1_end,
            structure::BasepairOP const & m2_end) {
        _link_motifs(m1, m1_end, m2, m2_end);
    }

    secondary_structure::PoseOP
    secondary_structure();

public:

    inline
    structure::ResidueOP
    get_residue(
            util::Uuid const & uuid) {
        auto r = structure::ResidueOP(nullptr);
        for (auto const & kv : motifs_) {
            r = kv.second->get_residue(uuid);
            if (r != nullptr) { return r; }
        }
        return nullptr;
    }

    inline
    structure::BasepairOP
    get_basepair(util::Uuid const & uuid) {
        if (all_bps_.find(uuid) != all_bps_.end()) {
            return all_bps_[uuid];
        } else {
            return nullptr;
        }
    }

    inline
    void
    to_pdb(
            String const & fname,
            int renumber = -1,
            int close_chains = 0,
            int conect_statements = 0) {
        return get_structure()->to_pdb(fname, renumber, close_chains, conect_statements);
    }

private:
    void
    _link_motifs(
            motif::MotifOP const &,
            structure::BasepairOP const &,
            motif::MotifOP const &,
            structure::BasepairOP const &);

    ChainNodes
    _get_end_nodes(
            ChainNodes const &,
            structure::BasepairOP const &);

    void
    _link_chains(
            ChainNodes &,
            ChainNodes &,
            MotifMergerType const &,
            MotifMergerType const &);

    void
    _connect_chains(
            ChainNode &,
            ChainNode &,
            int,
            int,
            MotifMergerType const &,
            MotifMergerType const &);

    void
    _build_structure();

    MotifMergerType
    _assign_merger_type(
            motif::MotifOP const &);

private:
    std::map<util::Uuid, structure::BasepairOP, util::UuidCompare> all_bps_;
    std::map<util::Uuid, motif::MotifOP, util::UuidCompare> motifs_;
    std::map<util::Uuid, util::Uuid, util::UuidCompare> res_overrides_, bp_overrides_;
    data_structure::graph::GraphStatic<ChainNodeData> graph_;
    int rebuild_structure_;
    structure::RNAStructureOP rna_structure_;

};

typedef std::shared_ptr<MotifMerger> MotifMergerOP;

}
#endif /* defined(__RNAMake__motif_merger__) */



































