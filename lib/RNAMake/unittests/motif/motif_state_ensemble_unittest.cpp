
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "math/numerical.h"
#include "motif/motif.h"
#include "motif/motif_state_ensemble.h"

TEST_CASE( "Test Motif State Ensembles", "[MotifStateEnsemble]" ) {
    auto path = base::motif_dirs() + "base.motif";
    auto m = motif::file_to_motif(path);
    
    auto mse = motif::MotifStateEnsemble(m->get_state());
    
    REQUIRE(mse.size() == 1);
    REQUIRE_THROWS_AS(mse.get_member(10), motif::MotifStateEnsembleException);
    REQUIRE_THROWS_AS(motif::MotifStateEnsemble(motif::MotifStateOPs(), Floats()), motif::MotifStateEnsembleException);
    REQUIRE_THROWS_AS(motif::MotifStateEnsemble(motif::MotifStateOPs{m->get_state()}, Floats()),
                      motif::MotifStateEnsembleException);

    
    SECTION("test copy constructor") {
        auto mse_copy = motif::MotifStateEnsemble(mse);
        
        REQUIRE(mse_copy.size() == 1);
    }

    SECTION("test stringifying motif state ensemble") {
        auto s = mse.to_str();
        auto mse_copy = motif::MotifStateEnsemble(s);
        
        REQUIRE(mse_copy.size() == 1);
    }
    
    
}