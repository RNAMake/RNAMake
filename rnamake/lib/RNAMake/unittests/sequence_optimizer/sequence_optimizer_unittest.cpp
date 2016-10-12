
//headers for testing
#include "../common.hpp"
#include "../tools/motif_graph_builder.hpp"


//RNAMake Headers
#include "base/settings.h"
#include "math/numerical.h"
#include "motif/motif.h"
#include "motif/motif_factory.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_topology.h"
#include "sequence_optimizer/sequence_optimizer_3d.hpp"

TEST_CASE( "Test Sequence Optimizer", "[SequenceOptimizer]" ) {
    
    SECTION("simplest test, optimize sequence to look like ideals in random structure") {
    
        auto builder = MotifGraphBuilder();
        auto mg = builder.build(2);
        auto m = RM::instance().motif("HAIRPIN.1C0A.0");
        mg->add_motif(m);
        mg->replace_ideal_helices();
    
        auto c = GraphtoTree();
        auto mt = c.convert(mg);
    
        auto so = SequenceOptimizer3D();
        so.set_option_value("verbose", false);
        auto sols = so.get_optimized_sequences(mt, mt->last_node()->data()->ends()[0],
                                               mt->last_node()->index(), 0);
    
        REQUIRE(sols.size() > 0);
    }
    
    
    SECTION("test optimizing miniTTR sequence") {
        
        auto base_path = base_dir() + "/rnamake/lib/RNAMake/apps/mini_ttr/resources/";
        RM::instance().add_motif(base_path+"GAAA_tetraloop");
        
        /*auto path = base_dir() + "/rnamake/unittests/test_problems/mini_ttr/sol.mg";
        auto lines = get_lines_from_file(path);
        auto mg = std::make_shared<MotifGraph>(lines[0], MotifGraphStringType::MG);
        mg->replace_ideal_helices();
        auto n1 = mg->get_node(1);
        auto n2 = n1->connections()[1]->partner(n1->index());
        
        auto c = GraphtoTree();
        auto mt = c.convert(mg, nullptr, -1, n2);
        mt->set_option_value("sterics", false);
        
        auto so = SequenceOptimizer3D();
        so.set_option_value("verbose", false);
        auto sols = so.get_optimized_sequences(mt, mt->last_node()->data()->ends()[0],
                                               mt->get_node(4)->index(), 1);*/

        //std::cout << sols.size() << std::endl;
        
        
    }
    
    
    
}