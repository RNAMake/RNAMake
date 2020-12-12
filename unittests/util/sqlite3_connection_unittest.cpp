

#include "../common.hpp"

#include "base/settings.h"
#include "util/sqlite3_connection.h"

TEST_CASE( "Test basic connection sqlite3 connection utilty" ) {
    
    SUBCASE("test catching nonexistant database files") {
        REQUIRE_THROWS_AS(util::Sqlite3Connection("test.db"), util::Sqlite3ConnectionException);
    }
    
    SUBCASE("require a database file to perform query") {
        auto sql_con = util::Sqlite3Connection();
        REQUIRE_THROWS_AS(sql_con.query("SELECT *"), util::Sqlite3ConnectionException);
    }
    
    auto path = base::resources_path()+"/motif_libraries_new/bp_steps.db";
    auto sql_con = util::Sqlite3Connection(path);

    SUBCASE("fetch first row of database") {
        auto row = sql_con.fetch_one("SELECT * from data_table");
        CHECK(row.size() == 5);
    }
    
    SUBCASE("count number of rows in database") {
        auto count = sql_con.count();
        CHECK(count > 0);
        
    }
    
}