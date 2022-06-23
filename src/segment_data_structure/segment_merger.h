////
//// Created by Joseph Yesselman on 1/23/18.
////
//
//#ifndef RNAMAKE_NEW_SEGMENT_MERGER_H
//#define RNAMAKE_NEW_SEGMENT_MERGER_H
//
//
//#include <resources/resource_manager.h>
//#include <secondary_structure/segment.h>
//#include <structure/segment.h>
//#include <segment_data_structure/segment_graph.h>
//
// namespace segment_data_structure {
//
// template <typename SegmentType>
// struct SegmentMergerResult {
//    inline
//    SegmentMergerResult(
//            std::shared_ptr<SegmentType> nsegment,
//            std::map<util::Uuid, util::Uuid> const & nres_uuid_map):
//            segment(nsegment),
//            res_uuid_map(nres_uuid_map) {}
//
//    std::shared_ptr<SegmentType> segment;
//    std::map<util::Uuid, util::Uuid> res_uuid_map;
//};
//
//
///*
// * Exception for segment merger
// */
// class SegmentMergerException : public std::runtime_error {
// public:
//    /**
//     * Standard constructor for SegmentMergerException
//     * @param   message   Error message for segment merger
//     */
//    SegmentMergerException(String const & message) :
//            std::runtime_error(message) {}
//};
//
// template<typename SegmentType, typename ChainType, typename ResType, typename
// BasepairType, typename AlignerType> class SegmentMerger { public:
//    struct ChainNodeData {
//        inline
//        ChainNodeData(
//                ChainType const & nchain,
//                util::Uuid const & nm_uuid,
//                bool nprime5_override,
//                bool nprime3_override):
//                chain(std::move(nchain)),
//                m_uuid(nm_uuid),
//                prime5_override(nprime5_override),
//                prime3_override(nprime3_override) {}
//
//        ChainType  chain;
//        util::Uuid m_uuid;
//        bool prime5_override;
//        bool prime3_override;
//    };
//
//    typedef data_structure::FixedEdgeUndirectedGraph<ChainNodeData>
//    ChainGraph;
//
// public:
//    enum class MergerType {
//        SPECIFIC_SEQUENCE, // sequence matters for these segments
//        NON_SPECIFIC_SEQUENCE // sequence does not matter such as idealized
//        helices
//    };
//
//    typedef std::shared_ptr<SegmentMergerResult<SegmentType> >
//    SegmentMergerResultOP;
//
// public:
//    SegmentMerger(
//            resources::ResourceManager const & rm):
//            rm_(rm),
//            end_chain_ids_1_(Indexes(2)),
//            end_chain_ids_2_(Indexes(2)),
//            res_uuid_map_(std::map<util::Uuid, util::Uuid>()) {}
//
//    ~SegmentMerger() {}
//
// public:
//    SegmentMergerResultOP
//    merge(
//            SegmentGraph<SegmentType, AlignerType> const & g,
//            String const & merged_name) {
//        res_uuid_map_ = std::map<util::Uuid, util::Uuid>();
//
//        for(auto const & n : g) {
//            if(!g.has_parent(n->index())) {
//                _add_segment(n->data());
//            }
//            else {
//                auto parent_index = g.get_parent_index(n->index());
//                auto parent_end_index = g.get_parent_end_index(n->index());
//                auto & segment_end = n->data().get_end(0);
//                auto & parent = g.get_segment(parent_index);
//                auto & parent_end = parent.get_end(parent_end_index);
//
//                _add_segment(n->data(), segment_end, parent, parent_end);
//            }
//        }
//
//        auto seg = _build_structure(g, merged_name);
//
//        return std::make_shared<SegmentMergerResult<SegmentType> >(seg,
//        res_uuid_map_);
//        //return SegmentMergerResultOP(nullptr);
//
//        //return std::make_shared<SegmentMergerResult>(seg, res_uuid_map_);
//
//    }
//
// protected:
//    void
//    _add_segment(
//            SegmentType const & sg) {
//        auto chains = sg.get_chains();
//        for(auto & c : *chains) {
//            chain_graph_.add_node(ChainNodeData(c, sg.get_uuid(), false,
//            false), 2);
//        }
//
//        for(auto it = sg.bps_begin(); it != sg.bps_end(); it++) {
//            all_bps_.push_back(&(*it));
//        }
//        chain_graph_.setup_transversal(0);
//
//    }
//
//    void
//    _add_segment(
//            SegmentType const & sg,
//            BasepairType const & sg_end,
//            SegmentType const & parent,
//            BasepairType const & parent_end) {
//        _add_segment(sg);
//        _link_segments(sg, sg_end, parent, parent_end);
//
//    }
//
//    void
//    _link_segments(
//            SegmentType const & sg,
//            BasepairType const & sg_end,
//            SegmentType const & parent,
//            BasepairType const & parent_end) {
//
//        _get_end_chains(sg_end, end_chain_ids_1_);
//        _get_end_chains(parent_end, end_chain_ids_2_);
//
//        auto sg_type_1 = _assign_merger_type(sg);
//        auto sg_type_2 = _assign_merger_type(parent);
//
//        if(sg_type_2 == MergerType::NON_SPECIFIC_SEQUENCE &&
//           sg_type_1 == MergerType::SPECIFIC_SEQUENCE) {
//            _link_chains(end_chain_ids_1_, end_chain_ids_2_, sg_type_1,
//            sg_type_2);
//        }
//        else {
//            _link_chains(end_chain_ids_2_, end_chain_ids_1_, sg_type_2,
//            sg_type_1);
//        }
//
//    }
//
//    void
//    _get_end_chains(
//            BasepairType const & end,
//            Indexes & end_chain_ids) {
//        end_chain_ids[0] = -1;
//        end_chain_ids[1] = -1;
//
//        for(auto const & n : chain_graph_) {
//            if     (n->data().chain.get_first().get_uuid() ==
//            end.get_res1_uuid() && end_chain_ids[0] == -1) {
//                end_chain_ids[0] = n->index();
//            }
//            else if(n->data().chain.get_first().get_uuid() ==
//            end.get_res1_uuid()) {
//                throw SegmentMergerException("end res1 is mapped to two
//                different chains!");
//            }
//
//            if     (n->data().chain.get_last().get_uuid() ==
//            end.get_res2_uuid() && end_chain_ids[1] == -1) {
//                end_chain_ids[1] = n->index();
//            }
//            else if(n->data().chain.get_last().get_uuid() ==
//            end.get_res2_uuid()) {
//                throw SegmentMergerException("end res2 is mapped to two
//                different chains!");
//            }
//        }
//
//        expects<SegmentMergerException>(
//                end_chain_ids[0] != -1 && end_chain_ids[1] != -1,
//                "did not find both chains for end residues");
//
//    }
//
//    MergerType
//    _assign_merger_type(
//            SegmentType const & sg) {
//        // not a helix probably cannot be replaced
//        if(sg.get_segment_type() != util::SegmentType::HELIX) {
//            return MergerType::SPECIFIC_SEQUENCE;
//        }
//
//        //ideal helices
//        if(sg.get_name_str().length() > 10 && sg.get_name_str().substr(0, 11)
//        == "HELIX.IDEAL") {
//            return MergerType::NON_SPECIFIC_SEQUENCE;
//        }
//        //flex helices
//        if(sg.get_name_str().length() > 9 && sg.get_name_str().substr(0, 10)
//        == "HELIX.FLEX") {
//            return MergerType::NON_SPECIFIC_SEQUENCE;
//        }
//
//        //bp step
//        if(sg.get_name_str().length() > 2 && sg.get_name_str().substr(0, 2) ==
//        "BP") {
//            return MergerType::NON_SPECIFIC_SEQUENCE;
//        }
//
//        return MergerType::SPECIFIC_SEQUENCE;
//    }
//
//    void
//    _link_chains(
//            Indexes dominant_indexes,
//            Indexes auxiliary_indexes,
//            MergerType sg_type_1,
//            MergerType sg_type_2) {
//
//        if(dominant_indexes[0] == dominant_indexes[1]) {
//            // dominant chain is a hairpin
//            _connect_chains(dominant_indexes[0], auxiliary_indexes[0], 1, 0,
//            sg_type_1, sg_type_2); _connect_chains(dominant_indexes[0],
//            auxiliary_indexes[1], 0, 1, sg_type_1, sg_type_2);
//        }
//        else if(auxiliary_indexes[0] == auxiliary_indexes[1]) {
//            // auxiliary chain is a hairpin
//            _connect_chains(dominant_indexes[1], auxiliary_indexes[0], 1, 0,
//            sg_type_1, sg_type_2); _connect_chains(dominant_indexes[0],
//            auxiliary_indexes[0], 0, 1, sg_type_1, sg_type_2);
//        }
//
//        else {
//            _connect_chains(dominant_indexes[1], auxiliary_indexes[0], 1, 0,
//            sg_type_1, sg_type_2); _connect_chains(dominant_indexes[0],
//            auxiliary_indexes[1], 0, 1, sg_type_1, sg_type_2);
//        }
//
//    }
//
//    void
//    _connect_chains(
//            Index d_index,
//            Index a_index,
//            Index d_end_index,
//            Index a_end_index,
//            MergerType sg_type_1,
//            MergerType sg_type_2) {
//
//        // auxiliary chain (a_index) first or last residue will be overritten
//        by the dominant chain
//        // during the final merger
//        if(a_end_index == 0) {
//            chain_graph_.get_node_data(a_index).prime5_override = 1;
//            res_uuid_map_[chain_graph_.get_node_data(a_index).chain.get_first().get_uuid()]
//            = \
//                          chain_graph_.get_node_data(d_index).chain.get_last().get_uuid();
//        }
//        else                 {
//            chain_graph_.get_node_data(a_index).prime3_override = 1;
//            res_uuid_map_[chain_graph_.get_node_data(a_index).chain.get_last().get_uuid()]
//            = \
//                          chain_graph_.get_node_data(d_index).chain.get_first().get_uuid();
//        }
//
//        chain_graph_.add_edge(data_structure::NodeIndexandEdge{d_index,
//        d_end_index},
//                              data_structure::NodeIndexandEdge{a_index,
//                              a_end_index});
//
//        // everything is okay no need for warnings
//        if(sg_type_1 == MergerType::NON_SPECIFIC_SEQUENCE ||
//           sg_type_2 == MergerType::NON_SPECIFIC_SEQUENCE) {
//            return;
//        }
//
//        if(a_end_index == 0 &&
//           chain_graph_.get_node_data(a_index).chain.get_first().get_name() !=
//           \ chain_graph_.get_node_data(d_index).chain.get_last().get_name())
//           {
//            _log_overriding_specific_squence_warning();
//        }
//        if(a_end_index == 1 &&
//           chain_graph_.get_node_data(d_index).chain.get_first().get_name() !=
//           \ chain_graph_.get_node_data(a_index).chain.get_last().get_name())
//           {
//            _log_overriding_specific_squence_warning();
//        }
//
//    }
//
//    void
//    _log_overriding_specific_squence_warning() {
//        LOGW << "Merging two chains with specific sequences this will remove
//        either the 5' or 3' of "; LOGW << "residue of a specified segments.
//        This likely is not what you wanted to do!!!";
//    }
//
//
//    void
//    _get_rna_residues_and_cutpoints(
//            std::vector<ResType> & res,
//            Cutpoints & cutpoints) {
//
//        auto start_indexes = Indexes();
//        for(auto const & n : chain_graph_) {
//            if(chain_graph_.edge_index_empty(n->index(), 0)) {
//            start_indexes.push_back(n->index()); }
//        }
//
//        for(auto const & si : start_indexes) {
//            auto * cur = &chain_graph_.get_node(si);
//            while(true) {
//                int start = 0;
//                int end = cur->data().chain.get_length()-1;
//                if(cur->data().prime5_override) { start = 1; }
//                if(cur->data().prime3_override) { end -= 1;  }
//                for(int i = start; i <= end; i++) {
//                    res.push_back(cur->data().chain.get_residue(i));
//                }
//
//                if(chain_graph_.edge_index_empty(cur->index(), 1)) { break; }
//                auto & edges = chain_graph_.get_node_edges(cur->index());
//                auto next_index = edges[1]->partner(cur->index());
//                cur = &chain_graph_.get_node(next_index);
//            }
//            cutpoints.push_back((int)res.size());
//        }
//
//    }
//
// public: // pure virtual needs to be different for each specialization
//
//    virtual
//    std::shared_ptr<SegmentType>
//    _build_structure(
//            SegmentGraph<SegmentType, AlignerType> const & g,
//            String const & merged_name) = 0;
//
//
// protected:
//    resources::ResourceManager const & rm_;
//    std::vector<BasepairType const *> all_bps_;
//    Indexes end_chain_ids_1_, end_chain_ids_2_;
//    std::map<util::Uuid, util::Uuid> res_uuid_map_;
//    ChainGraph chain_graph_;
//
//};
//
//}
//
//
// namespace secondary_structure {
//
// typedef segment_data_structure::SegmentMerger<Segment, Chain, Residue,
// Basepair, Aligner> _SegmentMerger;
//
// class SegmentMerger : public _SegmentMerger {
// public:
//    typedef _SegmentMerger BaseClass;
//
// public:
//    SegmentMerger(
//            resources::ResourceManager const & rm):
//            BaseClass(rm) {}
//
//
//    SegmentOP
//    _build_structure(
//            SegmentGraph const & g,
//            String const & merged_name) {
//
//        auto res = Residues();
//        auto chain_cuts = Cutpoints();
//        _get_rna_residues_and_cutpoints(res, chain_cuts);
//        auto res_uuids = std::map<util::Uuid, int>();
//        for(auto const & r : res) {
//            res_uuids[r.get_uuid()] = 1;
//        }
//
//        auto s = Structure(res, chain_cuts);
//
//        auto basepairs = Basepairs();
//
//        for(auto bp : all_bps_) {
//            if(res_uuids.find(bp->get_res1_uuid()) == res_uuids.end()) {
//            continue; } if(res_uuids.find(bp->get_res2_uuid()) ==
//            res_uuids.end()) { continue; } basepairs.push_back(*bp);
//        }
//
//        auto end_indexes = get_ends_from_basepairs(s, basepairs)->get_data();
//        auto end_ids = base::SimpleStringCOPs();
//        for(auto const & ei : end_indexes) {
//            auto end_id = generate_end_id(s, basepairs, basepairs[ei]);
//            end_ids.push_back(std::make_shared<base::SimpleString
//            const>(end_id));
//        }
//
//        auto name = std::make_shared<base::SimpleString const>(merged_name);
//        auto seg = std::make_shared<Segment>(s, basepairs, end_indexes,
//        end_ids, name,
//                                             util::SegmentType::SEGMENT, 0,
//                                             util::Uuid());
//        return seg;
//    }
//};
//
//
//}
//
//
// namespace structure {
//
// typedef segment_data_structure::SegmentMerger<Segment, Chain, Residue,
// Basepair, Aligner> _SegmentMerger;
//
// class SegmentMerger : public _SegmentMerger {
// public:
//    typedef _SegmentMerger BaseClass;
//
// public:
//    SegmentMerger(
//            resources::ResourceManager const & rm):
//            BaseClass(rm) {}
//
//    ~SegmentMerger() {}
//
// public:
//
//    SegmentOP
//    _build_structure(
//            SegmentGraph const & g,
//            String const & merged_name) {
//
//        auto res = Residues();
//        auto chain_cuts = Cutpoints();
//        _get_rna_residues_and_cutpoints(res, chain_cuts);
//        auto res_uuids = std::map<util::Uuid, int>();
//        for(auto const & r : res) {
//            res_uuids[r.get_uuid()] = 1;
//        }
//
//        auto s = Structure(res, chain_cuts);
//
//        auto protein_res = Residues();
//        auto protein_cutpoints = Cutpoints();
//
//        auto small_molecule_res = Residues();
//        auto small_molecule_cutpoints = Cutpoints();
//
//        for(auto const & n : g) {
//            for(auto it = n->data().protein_begin(); it !=
//            n->data().protein_end(); it += 1) {
//                protein_res.push_back(*it);
//                if (n->data().is_protein_residue_start_of_chain(*it)) {
//                    protein_cutpoints.push_back(protein_res.size());
//                }
//            }
//            if(protein_cutpoints.size() != 0 && protein_cutpoints.back() !=
//            protein_res.size()) {
//                protein_cutpoints.push_back(protein_res.size());
//            }
//
//            for(auto it = n->data().small_molecules_begin(); it !=
//            n->data().small_molecules_end(); it += 1) {
//                small_molecule_res.push_back(*it);
//                small_molecule_cutpoints.push_back(small_molecule_res.size());
//            }
//
//        }
//
//        auto basepairs = Basepairs();
//
//        for(auto bp : all_bps_) {
//            if(res_uuids.find(bp->get_res1_uuid()) == res_uuids.end()) {
//            continue; } if(res_uuids.find(bp->get_res2_uuid()) ==
//            res_uuids.end()) { continue; } basepairs.push_back(*bp);
//
//
//        }
//
//        return rm_.segment_from_components(merged_name, s, basepairs,
//                                           Structure(protein_res,
//                                           protein_cutpoints),
//                                           Structure(small_molecule_res,
//                                           small_molecule_cutpoints),
//                                           util::SegmentType::SEGMENT);
//
//    }
//
//
//};
//
//
//
//}
//
//#endif //RNAMAKE_NEW_SEGMENT_MERGER_H
