

//headers for testing
#include "../common.hpp"
#include "../tools/motif_graph_builder.hpp"

// Headers
#include "base/settings.h"
#include "math/numerical.h"
#include "motif/motif.h"
#include "motif/motif_factory.h"
#include "resources/resource_manager.h"
#include "motif_data_structure/motif_topology.h"
#include "sequence_optimization/sequence_optimizer_3d.hpp"

TEST_CASE( "Test Sequence Optimizer" ) {
    
    SUBCASE("test simple hairpin") {
        auto mg = std::make_shared<motif_data_structure::MotifGraph>();
        for(int i = 0; i < 10; i++) {
            mg->add_motif(resources::Manager::instance().motif("HELIX.IDEAL"));
        }
        mg->add_motif(resources::Manager::instance().motif("HAIRPIN.1GID.0"));
        auto end1 = std::make_shared<structure::Basepair>(*mg->get_node(9)->data()->ends()[1]);
        auto scorer = std::make_shared<sequence_optimization::ExternalTargetScorer>(
                mg->get_node(9)->data()->ends()[1]->state(), 8, 1, true);
        auto so = sequence_optimization::SequenceOptimizer3D();
        auto sols = so.get_optimized_sequences(mg, scorer);
        mg->replace_helical_sequence(sols[0]->sequence);
        auto end2 = mg->get_node(8)->data()->ends()[1];
        
        CHECK(mg->sequence() == sols[0]->sequence);
        auto diff = fabsf(sols[0]->dist_score - end1->diff(*end2));
        CHECK(diff < 0.1);
        
    }
    
    SUBCASE("test simple hairpin 2") {
        auto mg = std::make_shared<motif_data_structure::MotifGraph>();
        for(int i = 0; i < 10; i++) {
            mg->add_motif(resources::Manager::instance().motif("HELIX.IDEAL"));
        }
        mg->add_motif(resources::Manager::instance().motif("HAIRPIN.1GID.0"));
        auto end1 = std::make_shared<structure::Basepair>(*mg->get_node(9)->data()->ends()[1]);
        auto scorer = std::make_shared<sequence_optimization::ExternalTargetScorer>(
                mg->get_node(9)->data()->ends()[1]->state(), 8, 1, true);
        auto so = sequence_optimization::SequenceOptimizer3D();
        auto mg_opt = so.get_optimized_mg(mg, scorer);
        
        mg->replace_helical_sequence(mg_opt->sequence());
        CHECK(mg->sequence() == mg_opt->sequence());
        
        for(auto const & n : *mg) {
            CHECK(n->data()->name() == mg_opt->get_node(n->index())->data()->name());
        }
        
        auto atoms1 = mg->get_structure()->atoms();
        auto atoms2 = mg_opt->get_structure()->atoms();
        
        for(int i = 0; i < atoms1.size(); i++) {
            auto diff = atoms1[i]->coords().distance(atoms2[i]->coords());
            CHECK(diff < 0.01);
        }
        
    }
    
    SUBCASE("test optimizing miniTTR sequence") {
        auto path = base::base_dir() + "//unittests/test_problems/mini_ttr/sol.mg";
        auto lines =base::get_lines_from_file(path);
        auto mg = std::make_shared<motif_data_structure::MotifGraph>(lines[0],
                                                                     motif_data_structure::MotifGraphStringType::MG);
        
        mg->replace_ideal_helices();
        auto scorer = std::make_shared<sequence_optimization::InternalTargetScorer>(11, 1, 19, 1, false);
        auto so = sequence_optimization::SequenceOptimizer3D();
        auto mg_opt = so.get_optimized_mg(mg, scorer);
        
        mg->replace_helical_sequence(mg_opt->sequence());
        
        for(auto const & n : *mg) {
            CHECK(n->data()->name() == mg_opt->get_node(n->index())->data()->name());
        }
        
        auto atoms1 = mg->get_structure()->atoms();
        auto atoms2 = mg_opt->get_structure()->atoms();
        
        for(int i = 0; i < atoms1.size(); i++) {
            auto diff = atoms1[i]->coords().distance(atoms2[i]->coords());
            CHECK(diff < 0.01);
        }
        
    }


    // messed up backward compadility
    /*SUBCASE("test optimizing chip sequence") {
        auto base_path = base::base_dir() + "//unittests/resources/motif_graph";
        auto lines =base::get_lines_from_file(base_path+"/tecto_chip_only.mg");
        auto mg = std::make_shared<motif_data_structure::MotifGraph>(lines[0], motif_data_structure::MotifGraphStringType::MG);
        auto so = sequence_optimization::SequenceOptimizer3D();
        so.set_option_value("return_lowest", true);
        so.set_option_value("verbose", false);
        
        auto scorer = std::make_shared<InternalTargetScorer>(12, 1, 15, 1);
        auto mg_opt = so.get_optimized_mg(mg, scorer);
        
        mg->replace_helical_sequence(mg_opt->sequence());
        
        for(auto const & n : *mg) {
            CHECK(n->data()->name() == mg_opt->get_node(n->index())->data()->name());
        }
        
        std::cout << mg_opt->sequence() << std::endl;
        std::cout << mg_opt->dot_bracket() << std::endl;
    }*/
}
