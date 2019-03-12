
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif.h"
#include "resources/motif_sqlite_library.h"

TEST_CASE( "Test Motif Sqlite3 Library", "[MotifSqliteLibrary]" ) {
    
    SECTION("can load all the librarys listed in get_libnames") {
        for(auto const & kv : resources::MotifSqliteLibrary::get_libnames()) {
            REQUIRE_NOTHROW(resources::MotifSqliteLibrary(kv.first));
        }
    }
    
    SECTION("should return an error if not a valid library name") {
        REQUIRE_THROWS_AS(resources::MotifSqliteLibrary("fake"), resources::SqliteLibraryException);
        
    }
    
    SECTION("test ability to load all motifs from libaries") {
        for(auto const & kv : resources::MotifSqliteLibrary::get_libnames()) {
            auto mlib = resources::MotifSqliteLibrary(kv.first);
            REQUIRE_NOTHROW(mlib.load_all(10));
        }
    }
    
    SECTION("test individual queries") {
        auto mlib = resources::MotifSqliteLibrary("ideal_helices");
        REQUIRE(mlib.get("HELIX.IDEAL") != nullptr);
        REQUIRE(mlib.get("", "CC_LL_GG_RR") != nullptr);
        REQUIRE(mlib.get("", "", "A5-B7") != nullptr);
        REQUIRE(mlib.get("HELIX.IDEAL", "CC_LL_GG_RR") != nullptr);
        REQUIRE(mlib.get("HELIX.IDEAL", "", "A5-B7") != nullptr);
        REQUIRE(mlib.get("HELIX.IDEAL", "CC_LL_GG_RR", "A5-B7") != nullptr);
        REQUIRE(mlib.get("", "", "", "1") != nullptr);
        
        REQUIRE_THROWS_AS(mlib.get("TEST"), resources::SqliteLibraryException);
        REQUIRE_THROWS_AS(mlib.get("", "TEST"), resources::SqliteLibraryException);
        REQUIRE_THROWS_AS(mlib.get("", "", "TEST"), resources::SqliteLibraryException);
    }
    
    SECTION("test finding if a motif is contained in library") {
        auto mlib = resources::MotifSqliteLibrary("ideal_helices");
        
        REQUIRE(mlib.contains("HELIX.IDEAL") == 1);
        REQUIRE(mlib.contains("TEST") == 0);
        
    }

    // throws and error soemtimes but does not come up in the debugger?
    SECTION("test grabing more then one motif in the same query") {
        /*auto mlib = resources::MotifSqliteLibrary("twoway");
        auto m = mlib.get_random();
        auto motifs = mlib.get_multi(m->name());
         */
        
    }
    
}