//
// Created by Joseph Yesselman on 2/16/20.
//


#include "../common.hpp"

#include <base/settings.h>
#include <util/csv.h>

TEST_CASE( "Test csv parsing", "[CSV]" ) {

    SECTION("test reading a csv") {
        auto csv_reader = util::csv::Reader();
        auto csv_table = csv_reader.read_csv("test.csv");
        //std::cout << csv_table->num_rows() << std::endl;
        for(auto const & row : *csv_table) {
            std::cout << row.does_col_exist("test1") << std::endl;
        }
    }

    SECTION("test row topology") {
        auto row_top = util::csv::RowTopology(
                std::vector<DataType>{DataType::STRING, DataType::INT, DataType::FLOAT},
                {{"test1", 0}, {"test2", 1}, {"test3", 2}});

        REQUIRE(row_top.get_col_data_type("test1") == DataType::STRING);
        REQUIRE(row_top.get_col_data_type("test2") == DataType::INT);
        REQUIRE(row_top.get_col_data_type("test3") == DataType::FLOAT);

        REQUIRE_THROWS(row_top.get_col_data_type("fake"));

    }

    SECTION("test row") {
        auto row_top = util::csv::RowTopology(
                std::vector<DataType>{DataType::STRING, DataType::INT, DataType::FLOAT},
                {{"test1", 0}, {"test2", 1}, {"test3", 2}});

        auto data = std::vector<util::csv::Data>{
            {-1, -1, "I am a string"}, {1, -1, ""}, {-1, 44.2f, ""}
        };

        auto row = util::csv::Row(row_top, data);
        REQUIRE(row.get_string_val("test1") == "I am a string");
        REQUIRE(row.get_int_val("test2") == 1);
        REQUIRE(row.get_float_val("test3") == 44.2f);

        REQUIRE_THROWS(row.get_float_val("test1"));
        REQUIRE_THROWS(row.get_int_val("test3"));
    }

    /*SECTION("test catching nonexistant database files") {
        REQUIRE_THROWS_AS(util::Sqlite3Connection("test.db"), util::Sqlite3ConnectionException);
    }

    SECTION("require a database file to perform query") {
        auto sql_con = util::Sqlite3Connection();
        REQUIRE_THROWS_AS(sql_con.query("SELECT *"), util::Sqlite3ConnectionException);
    }

    auto path = base::resources_path()+"/motif_libraries_new/bp_steps.db";
    auto sql_con = util::Sqlite3Connection(path);

    SECTION("fetch first row of database") {
        auto row = sql_con.fetch_one("SELECT * from data_table");
        REQUIRE(row.size() == 5);
    }

    SECTION("count number of rows in database") {
        auto count = sql_con.count();
        REQUIRE(count > 0);

    }*/

}