////
//// Created by Joseph Yesselman on 5/31/18.
////
//
//#ifndef RNAMAKE_NEW_SEQUENCE_UPDATER_H
//#define RNAMAKE_NEW_SEQUENCE_UPDATER_H
//
//#include <resources/resource_manager.h>
//#include <segment_data_structure/segment_merger.h>
//
// namespace secondary_structure {
//
// class SequenceUpdater {
// public:
//    SequenceUpdater(
//            resources::ResourceManager const & rm):
//            sm_(SegmentMerger(rm)) {}
//
//    ~SequenceUpdater() {}
//
// public:
//    SegmentGraphOP
//    get_updated_graph(
//            SegmentGraph const & sg,
//            String const & seq) {
//        auto smr = sm_.merge(sg, "merged_graph");
//        merged_seg_ = smr->segment;
//        merged_seg_->set_sequence(seq);
//        auto new_sg = std::make_shared<SegmentGraph>(sg);
//        _update_graph_from_merged_segment(*new_sg, merged_seg_,
//        smr->res_uuid_map);
//
//        return new_sg;
//    }
//
//    void
//    update_graph(
//            SegmentGraph & sg,
//            String const & seq) {
//        auto smr = sm_.merge(sg, "merged_graph");
//        merged_seg_ = smr->segment;
//        merged_seg_->set_sequence(seq);
//        _update_graph_from_merged_segment(sg, merged_seg_, smr->res_uuid_map);
//
//    }
//
// private:
//    void
//    _update_graph_from_merged_segment(
//            SegmentGraph & sg,
//            SegmentOP seg,
//            std::map<util::Uuid, util::Uuid> const & res_uuid_map) {
//
//        auto c_uuid = util::Uuid();
//        for(auto & n : sg) {
//            int i = 0;
//            for(auto & r : n->data()) {
//                c_uuid = r.get_uuid();
//                if(res_uuid_map.find(c_uuid) != res_uuid_map.end()) {
//                    c_uuid = res_uuid_map.find(c_uuid)->second;
//            }
//                auto & m_r = merged_seg_->get_residue(c_uuid);
//                n->data().set_residue_identity(i, m_r.get_name());
//                i++;
//            }
//        }
//    }
//
//
// private:
//    SegmentMerger sm_;
//    SegmentOP merged_seg_;
//};
//
//
//}
//
// namespace structure {
//
// class SequenceUpdater {
// public:
//    SequenceUpdater(
//            resources::ResourceManager const & rm):
//            rm_(rm),
//            ss_su_(secondary_structure::SequenceUpdater(rm)) {}
//
//    ~SequenceUpdater() {}
//
// public:
//
//    SegmentGraphOP
//    get_updated_graph(
//            SegmentGraph const & sg,
//            String const & seq) {
//
//        auto ss_sg = get_secondary_structure_graph(sg);
//        auto sg_new = std::make_shared<SegmentGraph>(sg);
//        ss_su_.update_graph(*ss_sg, seq);
//        for(auto const & n : *ss_sg) {
//            //std::cout << n->data().get_sequence() << std::endl;
//            auto seg = rm_.get_segment(StringStringMap{{"end_id",
//            n->data().get_end_id(0)->get_str()}});
//            sg_new->replace_segment(n->index(), *seg, false);
//        }
//        return sg_new;
//    }
//
//
// private:
//    resources::ResourceManager const & rm_;
//    secondary_structure::SequenceUpdater ss_su_;
//
//};
//
//}
//
//
//
//#endif //RNAMAKE_NEW_SEQUENCE_UPDATER_H
