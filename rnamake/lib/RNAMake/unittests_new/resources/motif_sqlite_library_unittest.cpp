
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif.h"
#include "resources/motif_sqlite_library.h"

TEST_CASE( "Test Motif Sqlite3 Library", "[MotifSqliteLibrary]" ) {
    
    SECTION("can load all the librarys listed in get_libnames") {
        for(auto const & kv : MotifSqliteLibrary::get_libnames()) {
            REQUIRE_NOTHROW(MotifSqliteLibrary(kv.first));
        }
    }
    
    SECTION("should return an error if not a valid library name") {
        REQUIRE_THROWS_AS(MotifSqliteLibrary("fake"), SqliteLibraryException);
        
    }
    
    SECTION("test ability to load all motifs from libaries") {
        for(auto const & kv : MotifSqliteLibrary::get_libnames()) {
            auto mlib = MotifSqliteLibrary(kv.first);
            REQUIRE_NOTHROW(mlib.load_all(10));
        }
    }
    
    SECTION("test individual queries") {
        auto mlib = MotifSqliteLibrary("ideal_helices");
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
    
    SECTION("test finding if a motif is contained in library") {
        auto mlib = MotifSqliteLibrary("ideal_helices");
        
        REQUIRE(mlib.contains("HELIX.IDEAL") == 1);
        REQUIRE(mlib.contains("TEST") == 0);
        
    }

    // throws and error soemtimes but does not come up in the debugger?
    SECTION("test grabing more then one motif in the same query") {
        /*auto mlib = MotifSqliteLibrary("twoway");
        auto m = mlib.get_random();
        auto motifs = mlib.get_multi(m->name());
         */
        
    }
    
}