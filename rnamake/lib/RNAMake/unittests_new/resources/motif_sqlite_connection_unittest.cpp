
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif.h"
#include "resources/motif_sqlite_connection.h"

TEST_CASE( "Test Motif Sqlite3 Connection", "[MotifSqliteConnection]" ) {
    auto path = resources_path() + "motif_libraries_new/ideal_helices.db";
    auto conn = MotifSqliteConnection(path);
    
    SECTION("try querying for a motif") {
    
        auto query_str = String("SELECT * from data_table WHERE name='HELIX.IDEAL.3'");
    
        conn.query(query_str);
        auto row = conn.next();
    
        REQUIRE(row->name == "HELIX.IDEAL.3");
        conn.clear();
    }
    

    SECTION("invalid path") {
        REQUIRE_THROWS(MotifSqliteConnection("fake.db"));
    }
    
}