
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "math/numerical.h"
#include "motif/motif.h"
#include "motif/motif_ensemble.h"

TEST_CASE( "Test Motif Ensembles", "[MotifEnsemble]" ) {
    auto path = motif_dirs() + "base.motif";
    auto m = file_to_motif(path);

    auto me = MotifEnsemble(MotifOPs{m}, Floats{1});
    REQUIRE(me.size() == 1);
    
    SECTION("test copy constructor") {
        auto me_copy = MotifEnsemble(me);

        REQUIRE(me_copy.size() == 1);
        
    }
    
    SECTION("test stringifying motif ensemble") {
        auto s = me.to_str();
        auto rts = ResidueTypeSet();
        auto me_copy = MotifEnsemble(s, rts);
        
        REQUIRE(me_copy.size() == 1);

        
    }

}
