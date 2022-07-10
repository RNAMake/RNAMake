

#include "../common.hpp"

#include <base/paths.hpp>
#include <util/exception.hpp>
#include <util/sqlite/connection.hpp>

using namespace util::sqlite;

TEST_CASE("Test sqlite3 interface ") {
  SUBCASE("test database") {
    Database db;
    db.open(":memory:");
    CHECK(db.is_created());
    CHECK(db.is_open());
  }
  SUBCASE("test fields") {
    Field f("test", 1);
    CHECK(f.get_int() == 1);
    CHECK_THROWS_AS(double d = f.get_double(), util::SqliteException);
    CHECK_THROWS_AS(Blob b = f.get_blob(), util::SqliteException);
  }
  SUBCASE("test basic features of connection") {
    auto q = util::sqlite::Connection(":memory:");

    q.exec("CREATE TABLE data_table (word TEXT, id INT, PRIMARY KEY(id));");
    q.start_transaction();

    q.exec("INSERT INTO data_table( word, id) VALUES ('the_word','0');");
    q.exec("INSERT INTO data_table( word, id) VALUES ('the','1');");
    q.exec("INSERT INTO data_table( word, id) VALUES ('word','2');");

    q.commit_transaction();

    SUBCASE("test getting first row") {
      q.setup_row_iteration("SELECT * FROM data_table");
      auto &row = q.next();
      CHECK(row.size() == 2);
      CHECK(row[0].get_name() == "word");
      CHECK(row[0].get_type_str() == "the_word");

      int id = row[1].get_int();
      CHECK(id == 0);
    }

    SUBCASE("test going through all rows") {
      q.setup_row_iteration("SELECT * FROM data_table");
      int count = 0;
      while (q.has_more_rows()) {
        auto &row = q.next();
        count += 1;
      }
      CHECK(count == 3);
    }

    SUBCASE("test get table details") {
      auto tb = q.get_table_details("data_table");
      CHECK(tb->name() == "data_table");
      CHECK(tb->size() == 2);
      CHECK(tb->has_primary_key() == true);
      CHECK(tb->get_column(0).name == "word");
    }
  }
  SUBCASE("test create table") {
    auto conn = Connection(":memory:");
    auto td = TableDetails("data_table");
    td.add_column("word", "TEXT");
    td.add_column("id", "INT", true);
    util::sqlite::create_table(conn, td);

    SUBCASE("test reading table details from sqlite database") {
      auto new_td = *conn.get_table_details("data_table");
      CHECK(new_td.size() == 2);
      CHECK(new_td[0].name == "word");
      CHECK(new_td[0].type == "TEXT");
      CHECK(new_td[1].is_primary == true);
    }

    SUBCASE("test inserting multiple rows") {
      auto data =
          std::vector<Strings>{{"the_word", "0"}, {"the", "1"}, {"hello", "2"}};
      util::sqlite::insert_many(conn, "data_table", data);

      auto &row = conn.get_first_row("SELECT * FROM data_table");
      CHECK(row[0].get_name() == "word");
      CHECK(row[0].get_str() == "the_word");
    }
  }
  SUBCASE("test actual sql table") {
    String path =
        base::path::resources_path() + "/motif_libraries_new/ideal_helices.db";
    Connection q(path);
    auto const &row = q.get_first_row("SELECT * FROM data_table");
    CHECK(row[0].get_name() == "data");
    // example parsing existing motif
    Strings spl = base::string::split(row[0].get_str(), "&");
    CHECK(spl.size() == 11);
  }
  SUBCASE("test twoways") {
   // String path =
    //    base::path::resources_path() + "/motif_libraries_new/twoway.db";
    //Database db(path);
    //Connection q(db);
    /*try {
      auto const &row =
          q.setup_row_iteration("SELECT id ,data  FROM data_table WHERE "
                                "name='HELIX.IDEAL.2'");
    }
    catch(...) { }
    auto const &row =
        q.setup_row_iteration("SELECT id ,data  FROM data_table WHERE "
                              "name='HELIX.IDEAL.2'");  */
  }
}
