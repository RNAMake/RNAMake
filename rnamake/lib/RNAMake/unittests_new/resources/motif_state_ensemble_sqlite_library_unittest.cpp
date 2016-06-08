
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "resources/motif_state_ensemble_sqlite_library.h"

TEST_CASE( "Test Motif State Ensemble Library", "[MotifStateEnsembleSqliteLibrary]" ) {

    
    SECTION("can load all the librarys listed in get_libnames") {
        for(auto const & kv : MotifStateEnsembleSqliteLibrary::get_libnames()) {
            REQUIRE_NOTHROW(MotifStateEnsembleSqliteLibrary(kv.first));
            
        }
    }

    SECTION("should return an error if not a valid library name") {
        REQUIRE_THROWS_AS(MotifStateEnsembleSqliteLibrary("fake"), SqliteLibraryException);
        
    }
    
    SECTION("test ability to load all motifs states from libaries") {
        for(auto const & kv : MotifStateEnsembleSqliteLibrary::get_libnames()) {
            auto mlib = MotifStateEnsembleSqliteLibrary(kv.first);
            REQUIRE_NOTHROW(mlib.load_all(10));
        }
    }
    
    SECTION("test individual queries") {
        auto mlib = MotifStateEnsembleSqliteLibrary("bp_steps");
        REQUIRE(mlib.get("CC_LL_GG_RR") != nullptr);
        
        REQUIRE_THROWS_AS(mlib.get("TEST"), SqliteLibraryException);
        REQUIRE_THROWS_AS(mlib.get("", "TEST"), SqliteLibraryException);
    }


}