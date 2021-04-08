

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "resources/motif_ensemble_sqlite_connection.h"

TEST_CASE( "Test Motif Ensemble Sqlite Connection" ) {
    auto path = base::resources_path() + "motif_ensemble_libraries/bp_steps.db";
    auto conn = resources::MotifEnsembleSqliteConnection(path);
    
    SUBCASE("try querying for an ensemble") {
        auto query_str = String("SELECT * from data_table WHERE name='GG_LL_CC_RR'");
        
        conn.query(query_str);
        auto row = conn.next();
            
        CHECK(row->name == "GG_LL_CC_RR");
    }

    SUBCASE("invalid path") {
        REQUIRE_THROWS(resources::MotifEnsembleSqliteConnection("fake.db"));
    }
}