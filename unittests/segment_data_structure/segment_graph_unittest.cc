//
// Created by Joseph Yesselman on 1/22/18.
//


#include <iostream>
#include "../common.hpp"

#include <base/settings.h>
#include <base/log.h>
#include <resources/resource_manager.h>
#include <resources/segment_sqlite_library.h>
#include <segment_data_structure/segment_graph.h>

//TODO Why test Sqlite connections here?

TEST_CASE("Test Graph Data Structure") {

    SUBCASE("test ability to load all motifs from libaries") {
        for(auto const & kv : resources::SegmentSqliteLibrary::get_libnames()) {
            auto mlib = resources::SegmentSqliteLibrary(kv.first);
                    REQUIRE_NOTHROW(mlib.load_all(10));
        }
    }

    SUBCASE("basic tests") {
        auto rts = structure::ResidueTypeSet();
        auto db_path = base::resources_path() + "/motif_libraries/ideal_helices.db";
        auto seg_lib = resources::SegmentSqliteLibrary("data_table");
        auto seg1 = seg_lib.get();


        auto sg = structure::SegmentGraph();
        sg.add_segment(*seg1);


        for (int i = 0; i < 10; i++) {
            auto seg = seg_lib.get();
            sg.add_segment(*seg, i, sg.get_segment_end_name(i, 1));
        }

        CHECK(sg.get_num_segments() == 11);

        auto sg2 = structure::SegmentGraph(sg);
        CHECK(sg2.get_num_segments() == 11);

        sg.remove_segment(5);
        CHECK(sg2.get_num_segments() == 11);

        auto path = Indexes();
        auto target = Indexes{0, 1, 2, 3, 4, 6, 7, 8, 9, 10};
        for (auto const & n : sg) {
            path.push_back(n->index());
        }

        CHECK(path == target);

        // is copying messing up transveral
        path = Indexes();
        auto target_2 = Indexes{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        for (auto const & n : sg2) {
            path.push_back(n->index());
        }

        CHECK(path == target_2);

    }

//    SUBCASE("test add connection") {
//        auto sg = structure::SegmentGraph();
//        resources::ResourceManager rm;
//
//        auto seg1 = rm.get_segment(StringStringMap{{"name","HELIX.IDEAL.2"}});
//        sg.add_segment(*seg1);
//
//        for (int i = 0; i < 9; i++) {
//            auto seg = rm.get_segment(StringStringMap{{"name","HELIX.IDEAL.2"}});
//            sg.add_segment(*seg, i, sg.get_segment_end_name(i, 1));
//        }
//
//        sg.add_connection(data_structure::NodeIndexandEdge{0, 0},
//                          data_structure::NodeIndexandEdge{9, 1});
//
//        CHECK(sg.are_motifs_connected(0, 9) == true);
//
//    }

//    SUBCASE("test replace idealized helices") {
//        auto sg = structure::SegmentGraph();
//        resources::ResourceManager rm;
//        auto seg1 = rm.get_segment(StringStringMap{{"name","HELIX.IDEAL.10"}});
//        sg.add_segment(*seg1);
//
//        auto new_g = structure::convert_ideal_helices_to_basepair_steps(sg, rm);
//        CHECK(new_g->get_num_segments() == 11);
//    }
//
//    SUBCASE("test replace idealized helices with motifs") {
//        auto sg = structure::SegmentGraph();
//        resources::ResourceManager rm;
//
//        auto seg1 = rm.get_segment(StringStringMap{{"name","HELIX.IDEAL.5"}});
//        auto seg2 = rm.get_segment(StringStringMap{{"name","TWOWAY.1A34.0"}});
//        auto seg3 = rm.get_segment(StringStringMap{{"name","HELIX.IDEAL.5"}});
//        sg.add_segment(*seg1);
//        sg.add_segment(*seg2, 0, sg.get_segment_end_name(0, 1));
//        sg.add_segment(*seg3, 1, sg.get_segment_end_name(1, 1));
//        sg.add_connection(data_structure::NodeIndexandEdge{0, 0},
//                          data_structure::NodeIndexandEdge{2, 1});
//
//        auto new_g = structure::convert_ideal_helices_to_basepair_steps(sg, rm);
//        CHECK(new_g->get_num_segments() == 13);
//        CHECK(new_g->are_motifs_connected(0, 12));
//
//    }
//
//    SUBCASE("test replacing a segment") {
//        auto sg = structure::SegmentGraph();
//        auto sg2 = structure::SegmentGraph();
//        resources::ResourceManager rm;
//
//        auto seg1 = rm.get_segment(StringStringMap{{"name","HELIX.IDEAL.5"}});
//        auto seg2 = rm.get_segment(StringStringMap{{"name","TWOWAY.1A34.0"}});
//        auto seg3 = rm.get_segment(StringStringMap{{"name","HELIX.IDEAL.5"}});
//        auto seg4 = rm.get_segment(StringStringMap{{"name","HELIX.IDEAL.5"}});
//        sg.add_segment(*seg1);
//        sg.add_segment(*seg2, 0, sg.get_segment_end_name(0, 1));
//        sg.add_segment(*seg3, 1, sg.get_segment_end_name(1, 1));
//
//        sg.replace_segment(1, *seg4);
//
//        //compare to graph without replacement
//        sg2.add_segment(*seg1);
//        sg2.add_segment(*seg4, 0, sg2.get_segment_end_name(0, 1));
//        sg2.add_segment(*seg3, 1, sg2.get_segment_end_name(1, 1));
//
//        CHECK(sg.get_segment(2).get_end(1) == sg2.get_segment(2).get_end(1));
//    }
//
//    SUBCASE("test conversion to secondary structure") {
//        resources::ResourceManager rm;
//        auto seg1 = rm.get_segment(StringStringMap{{"name","HELIX.IDEAL.5"}});
//
//        auto ss = seg1->get_secondary_structure();
//        CHECK(ss->get_sequence() == seg1->get_sequence());
//        CHECK(ss->get_dot_bracket() == seg1->get_dot_bracket());
//
//        auto sg = structure::SegmentGraph();
//
//        auto seg2 = rm.get_segment(StringStringMap{{"name","TWOWAY.1A34.0"}});
//        auto seg3 = rm.get_segment(StringStringMap{{"name","HELIX.IDEAL.5"}});
//        sg.add_segment(*seg1);
//        sg.add_segment(*seg2, 0, sg.get_segment_end_name(0, 1));
//        sg.add_segment(*seg3, 1, sg.get_segment_end_name(1, 1));
//        sg.add_connection(data_structure::NodeIndexandEdge{0, 0},
//                          data_structure::NodeIndexandEdge{2, 1});
//
//        auto ss_sg = structure::get_secondary_structure_graph(sg);
//        CHECK(ss_sg->get_num_segments() == 3);
//        CHECK(ss_sg->are_motifs_connected(0, 2));
//
//    }

    //sg.write_nodes_to_pdbs("test");

}