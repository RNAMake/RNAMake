

#include "../common.hpp"

#include <base/paths.hpp>
#include <util/sqlite/connection.h>

using namespace util::sqlite;

TEST_CASE("Test sqlite3 interface ") {
  SUBCASE("test database") {
    Database db(":memory:");
    CHECK(db.is_created());
    CHECK(db.is_open());
  }

  /*SUBCASE("test basic features of connection") {
    auto db = util::sqlite::Database(":memory:");
    auto q = util::sqlite::Connection(db);

    q.exec("CREATE TABLE data_table (word TEXT, id INT, PRIMARY KEY(id));");
    q.start_transaction();

    q.exec("INSERT INTO data_table( word, id) VALUES ('the_word','0');");
    q.exec("INSERT INTO data_table( word, id) VALUES ('the','1');");
    q.exec("INSERT INTO data_table( word, id) VALUES ('word','2');");

    q.commit_transaction();

    SUBCASE("test getting first row") {
      q.setup_row_iteration("SELECT * FROM data_table");
      auto row = q.next();
      CHECK(row->at(0).get_name() == "word");
      CHECK(row->at(0).get_str() == "the_word");

      int id = row->at(1).get_int();
      CHECK(id == 0);

      //stop from weird casting ...
      //REQUIRE_THROWS_AS(int id = row->at(0), util::sqlite::SqliteException);
    }

    SUBCASE("test get table details") {
      auto tb = q.get_table_details("data_table");

    }
  }*/


}

/*TEST_CASE( "Test basic connection sqlite3 connection utilty" ) {
    
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
    
}      */