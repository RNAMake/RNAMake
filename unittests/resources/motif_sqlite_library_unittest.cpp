

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif.h"
#include "resources/motif_sqlite_library.h"

TEST_CASE( "Test Motif Sqlite3 Library" ) {
    
    SUBCASE("can load all the librarys listed in get_libnames") {
        for(auto const & kv : resources::MotifSqliteLibrary::get_libnames()) {
            REQUIRE_NOTHROW(resources::MotifSqliteLibrary(kv.first));
        }
    }
    
    SUBCASE("should return an error if not a valid library name") {
        REQUIRE_THROWS_AS(resources::MotifSqliteLibrary("fake"), resources::SqliteLibraryException);
        
    }
    
    SUBCASE("test ability to load all motifs from libaries") {
        for(auto const & kv : resources::MotifSqliteLibrary::get_libnames()) {
            auto mlib = resources::MotifSqliteLibrary(kv.first);
            REQUIRE_NOTHROW(mlib.load_all(10));
        }
    }
    
    SUBCASE("test individual queries") {
        auto mlib = resources::MotifSqliteLibrary("ideal_helices");
        CHECK(mlib.get("HELIX.IDEAL") != nullptr);
        CHECK(mlib.get("", "CC_LL_GG_RR") != nullptr);
        CHECK(mlib.get("", "", "A5-B7") != nullptr);
        CHECK(mlib.get("HELIX.IDEAL", "CC_LL_GG_RR") != nullptr);
        CHECK(mlib.get("HELIX.IDEAL", "", "A5-B7") != nullptr);
        CHECK(mlib.get("HELIX.IDEAL", "CC_LL_GG_RR", "A5-B7") != nullptr);
        CHECK(mlib.get("", "", "", "1") != nullptr);
        
        REQUIRE_THROWS_AS(mlib.get("TEST"), resources::SqliteLibraryException);
        REQUIRE_THROWS_AS(mlib.get("", "TEST"), resources::SqliteLibraryException);
        REQUIRE_THROWS_AS(mlib.get("", "", "TEST"), resources::SqliteLibraryException);
    }
    
    SUBCASE("test finding if a motif is contained in library") {
        auto mlib = resources::MotifSqliteLibrary("ideal_helices");
        
        CHECK(mlib.contains("HELIX.IDEAL") == 1);
        CHECK(mlib.contains("TEST") == 0);
        
    }

    // throws and error soemtimes but does not come up in the debugger?
    SUBCASE("test grabing more then one motif in the same query") {
        /*auto mlib = resources::MotifSqliteLibrary("twoway");
        auto m = mlib.get_random();
        auto motifs = mlib.get_multi(m->name());
         */
        
    }
    
}