
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "resources/resource_manager.h"
#include "resources/motif_sqlite_library.h"
#include "motif_data_structures/motif_merger.h"

TEST_CASE( "Test Mergering Motifs into single structure ", "[MotifMerger]" ) {
    auto m1 = RM::instance().motif("HELIX.IDEAL.3");
    auto m2 = RM::instance().motif("HELIX.IDEAL.3");
    m2->new_res_uuids();
    
    
    
}