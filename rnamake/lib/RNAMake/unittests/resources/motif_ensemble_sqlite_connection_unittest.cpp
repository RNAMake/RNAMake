
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "resources/motif_ensemble_sqlite_connection.h"

TEST_CASE( "Test Motif Ensemble Sqlite Connection", "[MotifEnsembleSqliteConnection]" ) {
    auto path = base::resources_path() + "motif_ensemble_libraries/bp_steps.db";
    auto conn = MotifEnsembleSqliteConnection(path);
    
    
    SECTION("try querying for an ensemble") {
        
        auto query_str = String("SELECT * from data_table WHERE name='GG_LL_CC_RR'");
        
        conn.query(query_str);
        auto row = conn.next();
            
        REQUIRE(row->name == "GG_LL_CC_RR");
    }
    
    
    SECTION("invalid path") {
        REQUIRE_THROWS(MotifEnsembleSqliteConnection("fake.db"));
    }


}