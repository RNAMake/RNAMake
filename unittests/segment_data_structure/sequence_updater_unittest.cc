//
// Created by Joseph Yesselman on 5/31/18.
//

//
// Created by Joseph Yesselman on 1/22/18.
//


#include <iostream>
#include "../common.hpp"

#include <base/settings.h>
#include <base/log.h>
#include <resources/resource_manager.h>
#include <segment_data_structure/sequence_updater.h>


TEST_CASE( "Test updating sequence in graphs ") {

    SUBCASE("test secondary structure sequence updating") {
        resources::ResourceManager rm;
        auto seg1 = rm.get_segment(StringStringMap{{"name", "HELIX.IDEAL"}});

        auto sg = structure::SegmentGraph();
        sg.add_segment(*seg1);

        for (int i = 0; i < 2; i++) {
            auto seg = rm.get_segment(StringStringMap{{"name", "HELIX.IDEAL"}});
            sg.add_segment(*seg, i, sg.get_segment_end_name(i, 1));
        }

        auto ss_sg = structure::get_secondary_structure_graph(sg);

        auto su = secondary_structure::SequenceUpdater(rm);
        auto ss_sg_new = su.get_updated_graph(*ss_sg, "AAAA&UUUU");
        for (auto const & n : *ss_sg_new) {
            CHECK(n->data().get_sequence() == "AA&UU");
            CHECK(n->data().get_end_id(0)->get_str() == "AA_LL_UU_RR");
        }
    }

    SUBCASE("test all atom sequence updating") {
        resources::ResourceManager rm;
        auto seg1 = rm.get_segment(StringStringMap{{"name", "HELIX.IDEAL"}});

        auto sg = structure::SegmentGraph();
        sg.add_segment(*seg1);

        for (int i = 0; i < 2; i++) {
            auto seg = rm.get_segment(StringStringMap{{"name", "HELIX.IDEAL"}});
            sg.add_segment(*seg, i, sg.get_segment_end_name(i, 1));
        }

        auto su = structure::SequenceUpdater(rm);
        auto sg_new = su.get_updated_graph(sg, String("AAAA&UUUU"));
        //sg_new->write_nodes_to_pdbs("test");
        auto end_id = sg_new->get_segment(0).get_end_id(0)->get_str();

        auto sg2 = structure::SegmentGraph();
        seg1 = rm.get_segment(StringStringMap{{"end_id", end_id}});
        sg2.add_segment(*seg1);
        for (int i = 0; i < 2; i++) {
            auto seg = rm.get_segment(StringStringMap{{"end_id", end_id}});
            sg2.add_segment(*seg, i, sg2.get_segment_end_name(i, 1));
        }

        CHECK(sg_new->get_segment(2).get_end(1).is_equal(sg2.get_segment(2).get_end(1), false));



    }

}
