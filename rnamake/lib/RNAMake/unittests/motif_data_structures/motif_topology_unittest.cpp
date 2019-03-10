
//headers for testing
#include "../common.hpp"
#include "../tools/motif_graph_builder.hpp"

//RNAMake Headers
#include "secondary_structure/util.h"
#include "motif_data_structures/motif_topology.h"
#include "resources/resource_manager.h"



TEST_CASE( "Test Changes in Motif Topology", "[MotifTopology]" ) {
    auto builder = MotifGraphBuilder();
    auto mg = builder.build(2);
    
    auto c = GraphtoTree();
    auto mt = c.convert(mg);
    REQUIRE(mt->size() == 4);
    
    mt = c.convert(mg, mg->last_node(), 1);
    REQUIRE(mt->size() == 4);
    
    mt = c.convert(mg, mg->get_node(0), 0);
    REQUIRE(mt->size() == 4);

    /*SECTION("start in the middle of the motif graph") {
    
        mt = c.convert(mg, mg->get_node(1), 1);
        REQUIRE(mt->size() == 2);
        
        mt = c.convert(mg, mg->get_node(1), 0);
        REQUIRE(mt->size() == 3);
        
    }
    
    mg->replace_ideal_helices();
    mt = c.convert(mg);
    
    REQUIRE(mt->size() == mg->size());

    auto build_points = mg->get_build_points();
    mt = c.convert(mg, build_points[0]->node, build_points[0]->end_index);
    REQUIRE(mt->size() == mg->size());
    
    auto dss = mg->designable_secondary_structure();
    fill_basepairs_in_ss(dss);
    mg->replace_helical_sequence(dss);
    
    mt = c.convert(mg);
    REQUIRE(mt->size() == mg->size());

    build_points = mg->get_build_points();
    mt = c.convert(mg, build_points[0]->node, build_points[0]->end_index);
    REQUIRE(mt->size() == mg->size());
    
    SECTION("test with minittr construct") {
        auto base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/mini_ttr/resources/";
        RM::instance().add_motif(base_path+"GAAA_tetraloop");
        
        auto path = base::base_dir() + "/rnamake/unittests/resources/motif_graph/mini_ttr.mg";
        auto lines =base::get_lines_from_file(path);
        auto mg = std::make_shared<MotifGraph>(lines[0], MotifGraphStringType::MG);

        auto c = GraphtoTree();
        auto mt = c.convert(mg);
        REQUIRE(mg->size() == mg->size());

        auto n = mg->get_node("GAAA_tetraloop");
        auto last_node = n->connections()[1]->partner(n->index());
        mt = c.convert(mg, nullptr, -1, last_node);
        REQUIRE(mg->size() == mg->size());

        //something strange happening here
        //mt = c.convert(mg, last_node);
        //REQUIRE(mg->size() == mg->size());
        
        REQUIRE(mt->connections().size() == 1);
     
        
        //std::cout << mt->to_pretty_str() << std::endl;
        
    }*/
    

    
}





















