//
// Created by Joseph Yesselman on 2019-03-30.
//

//
// Created by Joseph Yesselman on 3/23/19.
//



//headers for testing
#include "../common.hpp"
#include "../tools/motif_tree_builder.hpp"

//RNAMake Headers
#include "resources/resource_manager.h"

#include <base/log.h>
#include <motif_search/exhaustive/search.h>


TEST_CASE( "Test Searching Motif States" ) {
    auto & rm = resources::Manager::instance();
    auto start_ms = rm.motif_state("HELIX.IDEAL.2");

    SUBCASE("test enumerator") {
        using namespace motif_search::exhaustive;

        SUBCASE("test just flex helices") {
            auto flex_helices = resources::MotifStateSqliteLibrary("flex_helices");
            flex_helices.load_all();
            auto size = 0;
            for (auto const & ms : flex_helices) {
                if(ms->size() > 14) {
                    continue;
                }
                size += 1;
            }

            auto sol_template = motif_search::SolutionTopologyTemplate();
            sol_template.add_library("flex_helices");

            auto factory = motif_search::SolutionToplogyFactory();
            auto sol_toplogy = factory.generate_toplogy(sol_template);
            auto enumerator = MotifStateEnumerator(*sol_toplogy);
            auto count = 0;
            enumerator.start(start_ms->end_states()[1]);
            while (!enumerator.finished()) {
                enumerator.next();
                auto ms = enumerator.top_state();
                count += 1;
            }
            CHECK(size == count);

            // can repeat and not break?
            enumerator.start(start_ms->end_states()[1]);
            count = 0;
            while (!enumerator.finished()) {
                auto ms = enumerator.top_state();
                enumerator.next();
                count += 1;
            }
            CHECK(size == count);

        }

        SUBCASE("test flex helices and twoway") {
            auto flex_helices = resources::MotifStateSqliteLibrary("flex_helices");
            auto twoway = resources::MotifStateSqliteLibrary("twoway");
            flex_helices.load_all();
            twoway.load_all();
            auto size_1 = 0;
            auto size_2 = 0;
            for (auto const & ms : flex_helices) {
                if(ms->size() > 14) {
                    continue;
                }
                size_1 += 1;
            }
            for (auto const & ms : twoway) { size_2 += 1; }
            auto total_size = size_1 * size_2;

            auto sol_template = motif_search::SolutionTopologyTemplate();
            sol_template.add_library("flex_helices");
            sol_template.add_library("twoway", data_structure::NodeIndexandEdge{0, 1});

            auto factory = motif_search::SolutionToplogyFactory();
            auto sol_toplogy = factory.generate_toplogy(sol_template);
            auto enumerator = MotifStateEnumerator(*sol_toplogy);
            auto count = 0;
            enumerator.start(start_ms->end_states()[1]);
            while (!enumerator.finished()) {
                enumerator.next();
                count += 1;
            }

            CHECK(total_size == count);
            count = 0;
            enumerator.start(start_ms->end_states()[1]);
            while (!enumerator.finished()) {
                enumerator.next();
                count += 1;
                if (count == 100) { break; }
            }

            auto & motif_states = enumerator.all_states();

            auto aligner = motif::MotifStateAligner();
            auto ms1 = rm.motif_state(motif_states[0]->name(), "", motif_states[0]->end_names()[0]);
            auto ms2 = rm.motif_state(motif_states[1]->name(), "", motif_states[1]->end_names()[0]);

            aligner.get_aligned_motif_state(start_ms->end_states()[1], ms1);
            aligner.get_aligned_motif_state(ms1->end_states()[1], ms2);

            CHECK(ms1->end_states()[1]->diff(motif_states[0]->end_states()[1]) < 0.1);
            CHECK(ms2->end_states()[1]->diff(motif_states[1]->end_states()[1]) < 0.1);
        }
    }

    SUBCASE("test search") {
        SUBCASE("test flex helices only search") {
            using namespace motif_search::exhaustive;
            auto sol_template = motif_search::SolutionTopologyTemplate();
            sol_template.add_library("flex_helices");
            auto factory = motif_search::SolutionToplogyFactory();
            auto sol_toplogy = factory.generate_toplogy(sol_template);
            auto filter = std::make_shared<motif_search::NoExclusionFilter>();

            auto scorer = std::make_shared<DefaultScorer>();
            auto search = Search(scorer, *sol_toplogy, filter);

            auto start = start_ms->end_states()[0];
            auto end =  start_ms->end_states()[1];
            auto lookup = std::make_shared<util::StericLookupNew>();
            auto p = std::make_shared<motif_search::Problem>(start, end, lookup, true);

            search.setup(p);
            auto sol = search.next();
            //CHECK(sol != nullptr);
            //CHECK(sol->score < 10);

        }
    }
}




























