//
// Created by Joseph Yesselman on 3/23/19.
//


//headers for testing
#include "../common.hpp"
#include "../tools/motif_tree_builder.hpp"

//RNAMake Headers
#include "resources/resource_manager.h"

#include <base/log.h>
#include <motif_search/monte_carlo/search.h>

TEST_CASE( "Test Searching Motif States", "[MonteCarloSearch]" ) {
    auto & rm = resources::Manager::instance();

    SECTION("test scorer") {
        using namespace motif_search::monte_carlo;
        auto ms = rm.motif_state("HELIX.IDEAL.3");

        auto scorer = std::make_shared<GreedyScorer>();
        scorer->set_target(ms->end_states()[1], true);
        REQUIRE(scorer->score(*ms->end_states()[1]) == 0);

        scorer->set_target(ms->end_states()[1], false);
        REQUIRE(scorer->score(*ms->end_states()[1]) != 0);

    }

    SECTION("test solution topology template") {
        using namespace motif_search::monte_carlo;
        auto sol_template = motif_search::SolutionTopologyTemplate();
        sol_template.add_library("flex_helices");
        sol_template.add_library("twoway", data_structure::NodeIndexandEdge{0, 1});
        auto ms = rm.motif_state("HELIX.IDEAL.3");
        sol_template.add_motif_state(ms, data_structure::NodeIndexandEdge{1, 1});
        auto count = 0;
        for(auto const & n : sol_template) {
            count += 1;
        }
        REQUIRE(count == 3);
    }

    SECTION("test solution toplogy factory") {
        using namespace motif_search::monte_carlo;
        auto sol_template = motif_search::SolutionTopologyTemplate();
        sol_template.add_library("flex_helices");
        sol_template.add_library("twoway", data_structure::NodeIndexandEdge{0, 1});
        sol_template.add_library("flex_helices", data_structure::NodeIndexandEdge{1, 1});

        auto factory = motif_search::SolutionToplogyFactory();
        auto sol_toplogy = factory.generate_toplogy(sol_template);

        auto ms = rm.motif_state("HELIX.IDEAL.3");
        auto msg = sol_toplogy->initialize_solution(ms->end_states()[1]);

        auto diff = ms->end_states()[1]->diff(msg->get_node(1)->data()->get_end_state(0));
        REQUIRE(diff < 1);

        /*rm.get_motif_from_state(ms)->to_pdb("start.pdb");
        for(int i = 1; i < msg->size(); i++) {
            rm.get_motif_from_state(msg->get_node(i)->data()->cur_state)->to_pdb("test."+std::to_string(i)+".pdb");
        }*/
    }

    SECTION("test state") {
        using namespace motif_search::monte_carlo;
        auto mt = motif_data_structure::MotifTree();
        mt.add_motif(rm.motif("HELIX.IDEAL.3"));
        mt.add_motif(rm.motif("TWOWAY.2PN4.4"));
        mt.add_motif(rm.motif("HELIX.IDEAL.3"));

        auto start = mt.get_node(0)->data()->ends()[0]->state();
        auto end = mt.get_node(2)->data()->ends()[1]->state();
        auto lookup = std::make_shared<util::StericLookupNew>();
        auto p = std::make_shared<motif_search::Problem>(start, end, lookup, true);

        auto sol_template = motif_search::SolutionTopologyTemplate();
        sol_template.add_library("flex_helices");
        sol_template.add_library("twoway", data_structure::NodeIndexandEdge{0, 1});
        sol_template.add_library("flex_helices", data_structure::NodeIndexandEdge{1, 1});

        auto factory = motif_search::SolutionToplogyFactory();
        auto sol_toplogy = factory.generate_toplogy(sol_template);

        auto scorer = std::make_shared<GreedyScorer>();
        scorer->set_target(p->end, p->target_an_aligned_end);
        auto msg = sol_toplogy->initialize_solution(p->start);
        auto state = std::make_shared<State>(*msg, *sol_toplogy, p->lookup, scorer);

        REQUIRE(state->score == state->get_score());
        REQUIRE(state->steric_clash() == false);


        //auto search = Search(scorer, *sol_toplogy);
        //search.setup(p);

    }

}
































