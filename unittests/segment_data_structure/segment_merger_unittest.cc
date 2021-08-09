//
// Created by Joseph Yesselman on 1/22/18.
//


#include <iostream>
#include "../common.hpp"

#include <base/paths.h>
#include <base/log.h>
#include <resources/resource_manager.h>
#include <resources/segment_sqlite_library.h>
#include <segment_data_structure/segment_merger.h>


TEST_CASE( "Test Graph Data Structure ") {

    SUBCASE("test merger with structure") {

        resources::ResourceManager rm;
        auto seg1 = rm.get_segment(StringStringMap{{"name", "HELIX.IDEAL.2"}});

        auto sg = structure::SegmentGraph();
        sg.add_segment(*seg1);

        for (int i = 0; i < 10; i++) {
            auto seg = rm.get_segment(StringStringMap{{"name", "HELIX.IDEAL.2"}});
            sg.add_segment(*seg, i, sg.get_segment_end_name(i, 1));
        }

        auto sm = structure::SegmentMerger(rm);
        auto merged_results = sm.merge(sg, "merged_graph");
        auto merged_sg = merged_results->segment;

        CHECK(merged_results->res_uuid_map.size() == 20);
        CHECK(merged_sg->get_num_ends() == 2);
        CHECK(merged_sg->get_num_chains() == 2);
        CHECK(merged_sg->get_name_str() == "merged_graph");

        // copies normally
        auto merged_copy = structure::Segment(*merged_sg);
        CHECK(merged_copy.is_equal(*merged_sg));

        auto j = merged_sg->get_json();
        merged_copy = structure::Segment(j, rm.get_residue_type_set());
        CHECK(merged_copy.is_equal(*merged_sg, false));
    }

    SUBCASE("test merger on secondary structure") {
        resources::ResourceManager rm;
        auto seg1 = rm.get_segment(StringStringMap{{"name", "HELIX.IDEAL"}});

        auto sg = structure::SegmentGraph();
        sg.add_segment(*seg1);

        for (int i = 0; i < 2; i++) {
            auto seg = rm.get_segment(StringStringMap{{"name", "HELIX.IDEAL"}});
            sg.add_segment(*seg, i, sg.get_segment_end_name(i, 1));
        }

        auto ss_sg = structure::get_secondary_structure_graph(sg);
        auto sm = secondary_structure::SegmentMerger(rm);
        auto merged_results = sm.merge(*ss_sg, "merged_graph");
        auto merged_sg = merged_results->segment;

        CHECK(merged_sg->get_sequence() == "GGGG&CCCC");
        CHECK(merged_sg->get_dot_bracket() == "((((&))))");

    }


}






















