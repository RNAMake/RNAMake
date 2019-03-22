
//headers for testing
#include "../common.hpp"
#include "../tools/motif_tree_builder.hpp"

//RNAMake Headers
#include "resources/resource_manager.h"
#include "motif_search/motif_state_search.h"

#include <base/log.h>
#include <motif_search/path_finding/search.h>

TEST_CASE( "Test Searching Motif States", "[motif_search::MotifStateSearch]" ) {
    auto & rm = resources::Manager::instance();

    SECTION("test selector") {
        using namespace motif_search::path_finding;

        SECTION("test adding") {
            auto s = Selector();

            s.add("twoway");
            s.add(rm.motif("HELIX.IDEAL.3"));

            REQUIRE(s.size() == 2);

            REQUIRE_THROWS_AS(s.add("twoway"), SelectorException);
            REQUIRE_THROWS_AS(s.next(), SelectorException);

            s.start(0);
            REQUIRE(s.finished() == true);
        }

        SECTION("test connections") {
            auto s = Selector();

            s.add("twoway");
            s.add("ideal_helices");
            s.add(rm.motif("HELIX.IDEAL.3"));

            s.connect("twoway", "ideal_helices");

            SECTION("test iterating from twoway node") {

                s.start(0);
                REQUIRE(s.finished() == false);

                auto n = s.next();
                REQUIRE(n->name == "ideal_helices");
                REQUIRE(s.finished() == true);

            }

            SECTION("test iterating from ideal_helices node") {

                s.start(1);
                REQUIRE(s.finished() == false);

                auto n = s.next();
                REQUIRE(n->name == "twoway");
                REQUIRE(s.finished() == true);

            }

            s.connect("twoway", "HELIX.IDEAL.3");

            SECTION("test a second connection") {
                int count = 0;
                s.start(0);
                while(!s.finished()) {
                    count += 1;
                    s.next();
                }
                REQUIRE(count == 2);

            }

        }

        SECTION("test round robin selector") {
            auto s = RoundRobinSelector();

            s.add("ideal_helices");
            s.add("twoway");
            s.add("nway");
            s.add("flex_helices");

            REQUIRE(s.size() == 4);

            for(int i =0 ; i < 4; i++) {
                s.start(i);
                int count = 0;
                while (!s.finished()) {
                    count += 1;
                    s.next();
                }
                REQUIRE(count == 3);
            }
        }
    }

    SECTION("test nodes") {
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

        REQUIRE(n->score() == 0);
    }


    SECTION("test simple search") {
        using namespace motif_search::path_finding;
        base::init_logging(base::LogLevel::INFO);

        auto mt = motif_data_structure::MotifTree();
        mt.add_motif(rm.motif("HELIX.IDEAL.3"));

        auto start = mt.get_node(0)->data()->ends()[0]->state();
        auto end = mt.get_node(0)->data()->ends()[1]->state();
        auto lookup = std::make_shared<util::StericLookupNew>();
        auto p = std::make_shared<motif_search::Problem>(start, end, lookup, true);

        auto scorer = std::make_shared<GreedyBestFirst>();
        auto selector = default_selector();

        auto search = Search(scorer, selector);

        search.setup(p);
        auto sol = search.next();
        REQUIRE(sol != nullptr);

    }


    /*SECTION("test simple search") {
        auto mt = motif_data_structure::MotifTree();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.3");
        auto m2 = resources::Manager::instance().motif("TWOWAY.2PN4.4");
        auto m3 = resources::Manager::instance().motif("HELIX.IDEAL.3");
        mt.add_motif(m1);
        mt.add_motif(m2);
        mt.add_motif(m3);

        auto start = mt.get_node(0)->data()->ends()[0]->state();
        auto end = mt.get_node(2)->data()->ends()[1]->state();
        auto search = motif_search::MotifStateSearch();
        search.set_option_value("accept_score", 0.5f);
        search.set_option_value("max_node_level", 4);
        search.set_option_value("verbose", false);

        search.setup(start, end, true);
        auto sol = search.next();

        REQUIRE(sol != nullptr);
        REQUIRE(sol->score() < 0.1);

    }
    
    SECTION("test miniTTR") {
        resources::Manager::instance().add_motif(base::unittest_resource_dir() + "motif/GAAA_tetraloop");
        
        auto mt = motif_data_structure::MotifTree();
        auto m1 = resources::Manager::instance().motif("GAAA_tetraloop", "", "A229-A245");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.3");
        auto m3 = resources::Manager::instance().motif("HELIX.IDEAL.3");

        mt.add_motif(m1);
        mt.add_motif(m2);
        mt.add_motif(m3, 0);
        
        auto start = mt.get_node(1)->data()->ends()[1]->state();
        auto end = mt.get_node(2)->data()->ends()[1]->state();
        auto beads = mt.beads();
        auto centers = math::Points();

        for(auto const & b : beads) {
            if(b.btype() != structure::BeadType::PHOS) {
                centers.push_back(b.center());
            }
        }
        
        auto search = motif_search::MotifStateSearch();
        search.set_option_value("verbose", false);
        search.beads(centers);
        //search.set_option_value("max_node_level", 6);
        search.setup(start, end);
        auto sol = search.next();
        REQUIRE(sol->score() < 10);
        
        
    }*/
    
}
