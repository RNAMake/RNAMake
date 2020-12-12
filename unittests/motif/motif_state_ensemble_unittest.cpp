

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "math/numerical.h"
#include "motif/motif.h"
#include "motif/motif_state_ensemble.h"

#include <util/find_pair.h> // Is this import necessary?

TEST_CASE( "Test Motif State Ensembles" ) {
    auto path = base::motif_dirs() + "base.motif";
    auto m = motif::file_to_motif(path);
    
    auto mse = motif::MotifStateEnsemble(m->get_state());
    
    CHECK(mse.size() == 1);
    REQUIRE_THROWS_AS(mse.get_member(10), motif::MotifStateEnsembleException);
    REQUIRE_THROWS_AS(motif::MotifStateEnsemble(motif::MotifStateOPs(), Floats()), motif::MotifStateEnsembleException);
    REQUIRE_THROWS_AS(motif::MotifStateEnsemble(motif::MotifStateOPs{m->get_state()}, Floats()),
                      motif::MotifStateEnsembleException);

    
    SUBCASE("test copy constructor") {
        auto mse_copy = motif::MotifStateEnsemble(mse);
        
        CHECK(mse_copy.size() == 1);
    }

    SUBCASE("test stringifying motif state ensemble") {
        auto s = mse.to_str();
        auto mse_copy = motif::MotifStateEnsemble(s);
        
        CHECK(mse_copy.size() == 1);
    }
    
    
}