

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "resources/motif_state_sqlite_library.h"

TEST_CASE( "Test Motif State Sqlite3 Library" ) {
    
    SUBCASE("can load all the librarys listed in get_libnames") {
        for(auto const & kv : resources::MotifStateSqliteLibrary::get_libnames()) {
            REQUIRE_NOTHROW(resources::MotifStateSqliteLibrary(kv.first));
        }
    }
    
    SUBCASE("should return an error if not a valid library name") {
        REQUIRE_THROWS_AS(resources::MotifStateSqliteLibrary("fake"), resources::SqliteLibraryException);
        
    }
    
    SUBCASE("test ability to load all motifs states from libaries") {
        for(auto const & kv : resources::MotifStateSqliteLibrary::get_libnames()) {
            auto mlib = resources::MotifStateSqliteLibrary(kv.first);
            REQUIRE_NOTHROW(mlib.load_all(10));
        }
    }
    
    SUBCASE("test individual queries") {
        auto mlib = resources::MotifStateSqliteLibrary("ideal_helices");
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
    
    
    SUBCASE("test finding if a motif state is contained in library") {
        auto mlib = resources::MotifStateSqliteLibrary("ideal_helices");
        
        CHECK(mlib.contains("HELIX.IDEAL") == 1);
        CHECK(mlib.contains("TEST") == 0);
        
    }

    SUBCASE("test other libraries") {
        auto mlib = resources::MotifStateSqliteLibrary("nway");
        mlib.load_all(10);
        auto ct(0);  
        for(auto const & m : mlib) {
            //TODO add something else in here? I'm not sure what this is supposed to do  
            //std::cout << m->name() << std::endl;
            ++ct; 
        }
        //CHECK(ct == 10);
    }


}
