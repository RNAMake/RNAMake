

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif.h"
#include "resources/motif_sqlite_connection.h"

TEST_CASE( "Test Motif Sqlite3 Connection" ) {
    auto path = base::resources_path() + "motif_libraries_new/ideal_helices.db";
    auto conn = resources::MotifSqliteConnection(path);
    
    SUBCASE("try querying for a motif") {
        auto query_str = String("SELECT * from data_table WHERE name='HELIX.IDEAL.3'");
    
        conn.query(query_str);
        auto row = conn.next();
    
        CHECK(row->name == "HELIX.IDEAL.3");
        conn.clear();
    }

    SUBCASE("invalid path") {
        REQUIRE_THROWS(resources::MotifSqliteConnection("fake.db"));
    }
    
}