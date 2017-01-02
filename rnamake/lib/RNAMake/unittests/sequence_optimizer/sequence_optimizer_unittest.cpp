
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
    
    SECTION("test simple hairpin") {
        auto mg = std::make_shared<MotifGraph>();
        for(int i = 0; i < 10; i++) {
            mg->add_motif(RM::instance().motif("HELIX.IDEAL"));
        }
        mg->add_motif(RM::instance().motif("HAIRPIN.1GID.0"));
        auto end1 = std::make_shared<Basepair>(*mg->get_node(9)->data()->ends()[1]);
        auto scorer = std::make_shared<ExternalTargetScorer>(mg->get_node(9)->data()->ends()[1]->state(),
                                                             8, 1);
        auto so = SequenceOptimizer3D();
        auto sols = so.get_optimized_sequences(mg, scorer);
        mg->replace_helical_sequence(sols[0]->sequence);
        auto end2 = mg->get_node(8)->data()->ends()[1];
        
        REQUIRE(mg->sequence() == sols[0]->sequence);
        auto diff = fabsf(sols[0]->dist_score - end1->diff(*end2));
        REQUIRE(diff < 0.1);
        
    }
    
    SECTION("test simple hairpin 2") {
        auto mg = std::make_shared<MotifGraph>();
        for(int i = 0; i < 10; i++) {
            mg->add_motif(RM::instance().motif("HELIX.IDEAL"));
        }
        mg->add_motif(RM::instance().motif("HAIRPIN.1GID.0"));
        auto end1 = std::make_shared<Basepair>(*mg->get_node(9)->data()->ends()[1]);
        auto scorer = std::make_shared<ExternalTargetScorer>(mg->get_node(9)->data()->ends()[1]->state(),
                                                             8, 1);
        auto so = SequenceOptimizer3D();
        auto mg_opt = so.get_optimized_mg(mg, scorer);
        
        mg->replace_helical_sequence(mg_opt->sequence());
        REQUIRE(mg->sequence() == mg_opt->sequence());
        
        for(auto const & n : *mg) {
            REQUIRE(n->data()->name() == mg_opt->get_node(n->index())->data()->name());
        }
        
        auto atoms1 = mg->get_structure()->atoms();
        auto atoms2 = mg_opt->get_structure()->atoms();
        
        for(int i = 0; i < atoms1.size(); i++) {
            auto diff = atoms1[i]->coords().distance(atoms2[i]->coords());
            REQUIRE(diff < 0.01);
        }
        
    }
    
    SECTION("test optimizing miniTTR sequence") {
        auto path = base_dir() + "/rnamake/unittests/test_problems/mini_ttr/sol.mg";
        auto lines = get_lines_from_file(path);
        auto mg = std::make_shared<MotifGraph>(lines[0], MotifGraphStringType::MG);
        
        mg->replace_ideal_helices();
        auto scorer = std::make_shared<InternalTargetScorer>(11, 1, 19, 1);
        auto so = SequenceOptimizer3D();
        auto mg_opt = so.get_optimized_mg(mg, scorer);
        
        mg->replace_helical_sequence(mg_opt->sequence());
        
        for(auto const & n : *mg) {
            REQUIRE(n->data()->name() == mg_opt->get_node(n->index())->data()->name());
        }
        
        auto atoms1 = mg->get_structure()->atoms();
        auto atoms2 = mg_opt->get_structure()->atoms();
        
        for(int i = 0; i < atoms1.size(); i++) {
            auto diff = atoms1[i]->coords().distance(atoms2[i]->coords());
            REQUIRE(diff < 0.01);
        }
        
    }
    
    SECTION("test optimizing chip sequence") {
        auto base_path = base_dir() + "/rnamake/unittests/resources/motif_graph";
        auto lines = get_lines_from_file(base_path+"/tecto_chip_only.mg");
        auto mg = std::make_shared<MotifGraph>(lines[0], MotifGraphStringType::MG);
        mg->write_pdbs();
        auto so = SequenceOptimizer3D();
        so.set_option_value("return_lowest", true);
        so.set_option_value("verbose", false);
        
        auto scorer = std::make_shared<InternalTargetScorer>(12, 1, 15, 1);
        auto mg_opt = so.get_optimized_mg(mg, scorer);
        
        mg->replace_helical_sequence(mg_opt->sequence());
        
        for(auto const & n : *mg) {
            REQUIRE(n->data()->name() == mg_opt->get_node(n->index())->data()->name());
        }
        
        std::cout << mg_opt->sequence() << std::endl;
        std::cout << mg_opt->dot_bracket() << std::endl;
    }
}
