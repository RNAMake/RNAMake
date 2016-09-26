
#include "../common.hpp"

#include "base/settings.h"
#include "util/sqlite3_connection.h"

TEST_CASE( "Test basic connection sqlite3 connection utilty", "[Sqlite3Connection]" ) {
    
    SECTION("test catching nonexistant database files") {
        REQUIRE_THROWS_AS(Sqlite3Connection("test.db"), Sqlite3ConnectionException);
    }
    
    SECTION("require a database file to perform query") {
        auto sql_con = Sqlite3Connection();
        REQUIRE_THROWS_AS(sql_con.query("SELECT *"), Sqlite3ConnectionException);
    }
    
    auto path = resources_path()+"/motif_libraries_new/bp_steps.db";
    auto sql_con = Sqlite3Connection(path);

    SECTION("fetch first row of database") {
        auto row = sql_con.fetch_one("SELECT * from data_table");
        REQUIRE(row.size() == 5);
    }
    
    SECTION("count number of rows in database") {
        auto count = sql_con.count();
        REQUIRE(count > 0);
        
    }
    
}