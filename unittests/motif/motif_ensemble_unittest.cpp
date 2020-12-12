

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "math/numerical.h"
#include "motif/motif.h"
#include "motif/motif_ensemble.h"

TEST_CASE( "Test Motif Ensembles" ) {
    auto path = base::motif_dirs() + "base.motif";
    auto m = motif::file_to_motif(path);

    auto me = motif::MotifEnsemble(motif::MotifOPs{m}, Floats{1});
    CHECK(me.size() == 1);
    
    SUBCASE("test copy constructor") {
        auto me_copy = motif::MotifEnsemble(me);

        CHECK(me_copy.size() == 1);
        
    }
    
    SUBCASE("test stringifying motif ensemble") {
        auto s = me.to_str();
        auto rts = structure::ResidueTypeSet();
        auto me_copy = motif::MotifEnsemble(s, rts);
        
        CHECK(me_copy.size() == 1);

        
    }

}
