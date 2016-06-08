
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "resources/motif_state_sqlite_library.h"

TEST_CASE( "Test Motif State Sqlite3 Library", "[MotifStateSqliteLibrary]" ) {
    
    SECTION("can load all the librarys listed in get_libnames") {
        for(auto const & kv : MotifStateSqliteLibrary::get_libnames()) {
            REQUIRE_NOTHROW(MotifStateSqliteLibrary(kv.first));
    
        }
    }
    
    SECTION("should return an error if not a valid library name") {
        REQUIRE_THROWS_AS(MotifStateSqliteLibrary("fake"), SqliteLibraryException);
        
    }
    
    SECTION("test ability to load all motifs states from libaries") {
        for(auto const & kv : MotifStateSqliteLibrary::get_libnames()) {
            auto mlib = MotifStateSqliteLibrary(kv.first);
            REQUIRE_NOTHROW(mlib.load_all(10));
        }
    }
    
    SECTION("test individual queries") {
        auto mlib = MotifStateSqliteLibrary("ideal_helices");
        REQUIRE(mlib.get("HELIX.IDEAL") != nullptr);
        REQUIRE(mlib.get("", "CC_LL_GG_RR") != nullptr);
        REQUIRE(mlib.get("", "", "A5-B7") != nullptr);
        REQUIRE(mlib.get("HELIX.IDEAL", "CC_LL_GG_RR") != nullptr);
        REQUIRE(mlib.get("HELIX.IDEAL", "", "A5-B7") != nullptr);
        REQUIRE(mlib.get("HELIX.IDEAL", "CC_LL_GG_RR", "A5-B7") != nullptr);
        REQUIRE(mlib.get("", "", "", "1") != nullptr);
        
        REQUIRE_THROWS_AS(mlib.get("TEST"), SqliteLibraryException);
        REQUIRE_THROWS_AS(mlib.get("", "TEST"), SqliteLibraryException);
        REQUIRE_THROWS_AS(mlib.get("", "", "TEST"), SqliteLibraryException);
    }
    
    
    SECTION("test finding if a motif state is contained in library") {
        auto mlib = MotifStateSqliteLibrary("ideal_helices");
        
        REQUIRE(mlib.contains("HELIX.IDEAL") == 1);
        REQUIRE(mlib.contains("TEST") == 0);
        
    }


}