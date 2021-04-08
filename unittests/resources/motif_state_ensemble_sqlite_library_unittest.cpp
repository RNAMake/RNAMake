

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "resources/motif_state_ensemble_sqlite_library.h"

TEST_CASE( "Test Motif State Ensemble Library" ) {

    
    SUBCASE("can load all the librarys listed in get_libnames") {
        for(auto const & kv : resources::MotifStateEnsembleSqliteLibrary::get_libnames()) {
            REQUIRE_NOTHROW(resources::MotifStateEnsembleSqliteLibrary(kv.first));
        }
    }

    SUBCASE("should return an error if not a valid library name") {
        REQUIRE_THROWS_AS(resources::MotifStateEnsembleSqliteLibrary("fake"), resources::SqliteLibraryException);
        
    }
    
    SUBCASE("test ability to load all motifs states from libaries") {
        for(auto const & kv : resources::MotifStateEnsembleSqliteLibrary::get_libnames()) {
            auto mlib = resources::MotifStateEnsembleSqliteLibrary(kv.first);
            REQUIRE_NOTHROW(mlib.load_all(10));
        }
    }
    
    SUBCASE("test individual queries") {
        auto mlib = resources::MotifStateEnsembleSqliteLibrary("bp_steps");
        CHECK(mlib.get("CC_LL_GG_RR") != nullptr);
        
        REQUIRE_THROWS_AS(mlib.get("TEST"), resources::SqliteLibraryException);
        REQUIRE_THROWS_AS(mlib.get("", "TEST"), resources::SqliteLibraryException);
    }


}
