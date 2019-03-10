
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif_scorer.h"

TEST_CASE( "Test Secondary Structure scoring of motifs", "[MotifScorer]" ) {
    
    auto path = base::motif_dirs() + "base.motif";
    auto m = file_to_motif(path);
    path = base::motif_dirs() + "ref.motif";
    auto ref_m = file_to_motif(path);;

    auto scorer = MotifScorer();
    REQUIRE(scorer.score(ref_m) == -1.9f);
    REQUIRE(scorer.score(m) == -9.5f);

}