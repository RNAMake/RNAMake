//
// Created by Joseph Yesselman on 2/16/20.
//



#include "../common.hpp"

#include <base/paths.hpp>
#include <util/csv.h>

TEST_CASE( "Test csv parsing" ) {
  SUBCASE("test reading a csv") {
    auto path =  base::path::unittest_resource_path() + "/base/test.csv";
    auto csv_reader = io::CSVReader<3>(path);
    String test1, test2, test3;
    int count = 0;
    while(csv_reader.read_row(test1, test2, test3)) {
      count += 1;
    }
    CHECK(count == 3);
  }

    /*SUBCASE("test reading a csv") {
        auto csv_reader = util::csv::Reader();
        auto path =  base::unittest_resource_dir() + "/base/test.csv";
        auto csv_table = csv_reader.read_csv(path);
        CHECK(csv_table->num_rows() == 2);
        CHECK(csv_table->does_col_exist("test1"));
        CHECK(!csv_table->does_col_exist("fake"));

        auto col_1_data = Strings{
            "I am a string",
            "I am a string 2"
        };

        auto col_2_data = Ints{1, 2};
        auto col_3_data = Floats{44.2, 45.2};

        auto pos = 0;
        for(auto const & row : *csv_table) {
            CHECK(row.does_col_exist("test1"));
            CHECK(row.get_string_val("test1") == col_1_data[pos]);
            CHECK(row.get_int_val("test2") == col_2_data[pos]);
            CHECK(row.get_float_val("test3") == col_3_data[pos]);
            pos += 1;
        }
    }

    SUBCASE("error catching for table production") {
        auto csv_reader = util::csv::Reader();
        // not a real file
        auto fake_path = base::unittest_resource_dir() + "/base/fake.csv";
        REQUIRE_THROWS(csv_reader.read_csv(fake_path));

        // is a file but nothing in it
        auto path_1 = base::unittest_resource_dir() + "/base/no_lines.csv";
        REQUIRE_THROWS(csv_reader.read_csv(path_1));

        // there are col names but no data
        auto path_2 = base::unittest_resource_dir() + "/base/no_data_lines.csv";
        REQUIRE_THROWS(csv_reader.read_csv(path_2));

        // incorrect number of data cols vs names
        auto path_3 = base::unittest_resource_dir() + "/base/incorrect_num_of_data_lines.csv";
        REQUIRE_THROWS(csv_reader.read_csv(path_3));

        //
    }

    SUBCASE("test row topology") {
        auto data_types = DataTypes{DataType::STRING, DataType::INT, DataType::FLOAT};
        auto col_map = std::map<String, int>{{"test1", 0}, {"test2", 1}, {"test3", 2}};
        auto row_top = util::csv::RowTopology(data_types, col_map);

        CHECK(row_top.get_col_data_type("test1") == DataType::STRING);
        CHECK(row_top.get_col_data_type("test2") == DataType::INT);
        CHECK(row_top.get_col_data_type("test3") == DataType::FLOAT);

        REQUIRE_THROWS(row_top.get_col_data_type("fake"));

    }

    SUBCASE("test row") {
        auto data_types = DataTypes{DataType::STRING, DataType::INT, DataType::FLOAT};
        auto col_map = std::map<String, int>{{"test1", 0}, {"test2", 1}, {"test3", 2}};
        auto row_top = util::csv::RowTopology(data_types, col_map);

        auto data = std::vector<util::csv::Data>{
            {-1, -1, "I am a string"}, {1, -1, ""}, {-1, 44.2f, ""}
        };

        auto row = util::csv::Row(row_top, data);
        CHECK(row.get_string_val("test1") == "I am a string");
        CHECK(row.get_int_val("test2") == 1);
        CHECK(row.get_float_val("test3") == 44.2f);

        REQUIRE_THROWS(row.get_float_val("test1"));
        REQUIRE_THROWS(row.get_int_val("test3"));
    }
     */

}