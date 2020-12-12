

//headers for testing
#include "../common.hpp"
#include "../tools/motif_tree_builder.hpp"

//RNAMake Headers
#include "resources/resource_manager.h"

#include <base/log.h>
#include <motif_search/path_finding/search.h>

TEST_CASE( "Test Searching Motif States" ) {
   auto & rm = resources::Manager::instance();

    SUBCASE("test selector") {
        using namespace motif_search::path_finding;

        SUBCASE("test adding") {
            auto s = Selector();

            s.add("twoway");
            s.add(rm.motif("HELIX.IDEAL.3"));

            CHECK(s.size() == 2);

            REQUIRE_THROWS_AS(s.add("twoway"), SelectorException);
            REQUIRE_THROWS_AS(s.next(), SelectorException);

            s.start(0);
            CHECK(s.finished() == true);
        }

        SUBCASE("test connections") {
            auto s = Selector();

            s.add("twoway");
            s.add("ideal_helices");
            s.add(rm.motif("HELIX.IDEAL.3"));

            s.connect("twoway", "ideal_helices");

            SUBCASE("test iterating from twoway node") {

                s.start(0);
                CHECK(s.finished() == false);

                auto n = s.next();
                CHECK(n->name == "ideal_helices");
                CHECK(s.finished() == true);

            }

            SUBCASE("test iterating from ideal_helices node") {

                s.start(1);
                CHECK(s.finished() == false);

                auto n = s.next();
                CHECK(n->name == "twoway");
                CHECK(s.finished() == true);

            }

            s.connect("twoway", "HELIX.IDEAL.3");

            SUBCASE("test a second connection") {
                int count = 0;
                s.start(0);
                while(!s.finished()) {
                    count += 1;
                    s.next();
                }
                CHECK(count == 2);

            }

        }

        SUBCASE("test round robin selector") {
            auto s = RoundRobinSelector();

            s.add("ideal_helices");
            s.add("twoway");
            s.add("nway");
            s.add("flex_helices");

            CHECK(s.size() == 4);

            for(int i =0 ; i < 4; i++) {
                s.start(i);
                int count = 0;
                while (!s.finished()) {
                    count += 1;
                    s.next();
                }
                CHECK(count == 3);
            }
        }
    }

    SUBCASE("test nodes") {
        using namespace motif_search::path_finding;
        auto ms = rm.motif_state("HELIX.IDEAL.3");
        auto n1 = std::make_shared<Node>(ms, nullptr, 0, 0, 0, 0);
        auto n2 = std::make_shared<Node>(ms, n1, 10, 0, 0, 0);

        auto comparer = NodeCompare();
        comparer(n1, n2);

        auto queue = NodeQueue();
        queue.push(n1);
        queue.push(n2);

        auto n = queue.top();

        CHECK(n->score() == 0);
    }


    SUBCASE("test simple search") {
        using namespace motif_search::path_finding;
        //base::init_logging(base::LogLevel::INFO);

        auto mt = motif_data_structure::MotifTree();
        mt.add_motif(rm.motif("HELIX.IDEAL.3"));

        auto start = mt.get_node(0)->data()->ends()[0]->state();
        auto end = mt.get_node(0)->data()->ends()[1]->state();
        auto lookup = std::make_shared<util::StericLookupNew>();
        auto p = std::make_shared<motif_search::Problem>(start, end, lookup, true);

        auto scorer = std::make_shared<GreedyScorer>();
        auto selector = default_selector();
        auto filter = std::make_shared<motif_search::NoExclusionFilter>();

        auto search = Search(scorer, selector, filter);

        search.setup(p);
        auto sol = search.next();
        CHECK(sol != nullptr);

    }

    SUBCASE("test contrained searches") {
        using namespace motif_search::path_finding;
        auto mt = motif_data_structure::MotifTree();
        mt.add_motif(rm.motif("HELIX.IDEAL.3"));
        mt.add_motif(rm.motif("TWOWAY.2PN4.4"));
        mt.add_motif(rm.motif("HELIX.IDEAL.3"));

        auto start = mt.get_node(0)->data()->ends()[0]->state();
        auto end = mt.get_node(2)->data()->ends()[1]->state();
        auto lookup = std::make_shared<util::StericLookupNew>();
        auto p = std::make_shared<motif_search::Problem>(start, end, lookup, true);

        auto scorer = std::make_shared<GreedyScorer>();
        auto selector = default_selector();
        auto filter = std::make_shared<motif_search::NoExclusionFilter>();

        auto search = Search(scorer, selector, filter);
        search.setup(p);
        auto sol = search.next();
        CHECK(sol != nullptr);

        // max number of motifs
        search = Search(scorer, selector, filter);
        search.setup(p);
        search.set_option_value("max_node_level", 4);
        sol = search.next();
        CHECK(sol != nullptr);
        CHECK(sol->graph->size() <= 3);

        // max number of residues in solution
        search = Search(scorer, selector, filter);
        search.setup(p);
        search.set_option_value("max_size", 30);
        sol = search.next();
        CHECK(sol != nullptr);

        auto res_count = 2;
        for(auto const & n : *sol->graph) {
            res_count += n->data()->cur_state->size() - 2;
        }
        CHECK(res_count <= 30);

        // solution must end in helix
        search = Search(scorer, selector, filter);
        search.setup(p);
        search.set_option_value("helix_end", true);
        sol = search.next();
        CHECK(sol != nullptr);
        CHECK(sol->graph->last_node()->data()->cur_state->name()[0] == 'H');

        // should return no solution
        search = Search(scorer, selector, filter);
        search.setup(p);
        search.set_option_value("max_node_level", 2);
        sol = search.next();
        CHECK(sol == nullptr);

        // should return the best solution but not within cutoff
        search = Search(scorer, selector, filter);
        search.setup(p);
        search.set_option_value("max_node_level", 2);
        search.set_option_value("return_best", true);
        sol = search.next();
        CHECK(sol != nullptr);
        CHECK(sol->score > 10);

    }

    SUBCASE("test solution filtering") {
        auto filter = std::make_shared<motif_search::RemoveDuplicateHelices>();

        auto m_names_1 = Strings{"HELIX.FLEX.5.5", "TWOWAY.2PN4.4", "HELIX.FLEX.6.5"};
        auto m_names_2 = Strings{"HELIX.FLEX.5.2", "TWOWAY.2PN4.4", "HELIX.FLEX.6.2"};
        auto m_names_3 = Strings{"HELIX.FLEX.6.2", "TWOWAY.2PN4.4", "HELIX.FLEX.5.2"};

        CHECK(filter->accept(m_names_1) == true);
        // can't accept the same solution twice
        CHECK(filter->accept(m_names_1) == false);
        // can't accept solution with the same helices
        CHECK(filter->accept(m_names_2) == false);
        // should accept new solution
        CHECK(filter->accept(m_names_3) == true);

    }

    SUBCASE("test miniTTR") {
        using namespace motif_search::path_finding;
        base::init_logging(base::LogLevel::WARN);
        rm.add_motif(base::unittest_resource_dir() + "motif/GAAA_tetraloop");


        auto mt = motif_data_structure::MotifTree();
        auto m1 = rm.motif("GAAA_tetraloop", "", "A229-A245");
        auto m2 = rm.motif("HELIX.IDEAL.3");
        auto m3 = rm.motif("HELIX.IDEAL.3");

        mt.add_motif(m1);
        mt.add_motif(m2);
        mt.add_motif(m3, 0);

        auto start = mt.get_node(1)->data()->ends()[1]->state();
        auto end = mt.get_node(2)->data()->ends()[1]->state();
        auto lookup = std::make_shared<util::StericLookupNew>();
        auto p = std::make_shared<motif_search::Problem>(start, end, lookup, false);

        auto beads = mt.beads();
        auto centers = math::Points();

        for (auto const &b : beads) {
            if (b.btype() != structure::BeadType::PHOS) {
                centers.push_back(b.center());
            }
        }
        //lookup->add_points(centers);

        auto scorer = std::make_shared<AstarScorer>();
        auto selector = default_selector();
        auto filter = std::make_shared<motif_search::NoExclusionFilter>();

        auto search = Search(scorer, selector, filter);
        search.set_option_value("accept_score", 5);
        search.setup(p);
        auto sol = search.next();
        //CHECK(sol != nullptr);

    }

}
