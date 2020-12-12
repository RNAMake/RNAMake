

//headers for testing
#include "../common.hpp"
#include "../tools/motif_graph_builder.hpp"

//RNAMake Headers
#include "secondary_structure/util.h"
#include "motif_data_structure/motif_topology.h"
#include "resources/resource_manager.h"

#include <util/find_pair.h>

TEST_CASE( "Test Changes in Motif Topology" ) {
    auto builder = MotifGraphBuilder();
    auto mg = builder.build(2);
    
    auto c = motif_data_structure::GraphtoTree();
    auto mt = c.convert(mg);
    CHECK(mt->size() == 4);
    
    mt = c.convert(mg, mg->last_node(), 1);
    CHECK(mt->size() == 4);
    
    mt = c.convert(mg, mg->get_node(0), 0);
    CHECK(mt->size() == 4);

    /*SUBCASE("start in the middle of the motif graph") {
    
        mt = c.convert(mg, mg->get_node(1), 1);
        CHECK(mt->size() == 2);
        
        mt = c.convert(mg, mg->get_node(1), 0);
        CHECK(mt->size() == 3);
        
    }
    
    mg->replace_ideal_helices();
    mt = c.convert(mg);
    
    CHECK(mt->size() == mg->size());

    auto build_points = mg->get_build_points();
    mt = c.convert(mg, build_points[0]->node, build_points[0]->end_index);
    CHECK(mt->size() == mg->size());
    
    auto dss = mg->designable_secondary_structure();
    fill_basepairs_in_ss(dss);
    mg->replace_helical_sequence(dss);
    
    mt = c.convert(mg);
    CHECK(mt->size() == mg->size());

    build_points = mg->get_build_points();
    mt = c.convert(mg, build_points[0]->node, build_points[0]->end_index);
    CHECK(mt->size() == mg->size());
    
    SUBCASE("test with minittr construct") {
        auto base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/mini_ttr/resources/";
        resources::Manager::instance().add_motif(base_path+"GAAA_tetraloop");
        
        auto path = base::base_dir() + "/rnamake/unittests/resources/motif_graph/mini_ttr.mg";
        auto lines =base::get_lines_from_file(path);
        auto mg = std::make_shared<motif_data_structure::MotifGraph>(lines[0], motif_data_structure::MotifGraphStringType::MG);

        auto c = motif_data_structure::GraphtoTree();
        auto mt = c.convert(mg);
        CHECK(mg->size() == mg->size());

        auto n = mg->get_node("GAAA_tetraloop");
        auto last_node = n->connections()[1]->partner(n->index());
        mt = c.convert(mg, nullptr, -1, last_node);
        CHECK(mg->size() == mg->size());

        //something strange happening here
        //mt = c.convert(mg, last_node);
        //CHECK(mg->size() == mg->size());
        
        CHECK(mt->connections().size() == 1);
     
        
        //std::cout << mt->to_pretty_str() << std::endl;
        
    }*/
    

    
}





















