//
// Created by Joseph Yesselman on 2/16/20.
//

//
// Created by Joseph Yesselman on 3/10/19.
//


#include <iostream>
#include "../common.hpp"
#include <base/string.h>

TEST_CASE("test string functions", "[String]" ) {
    SECTION("test determine_string_data_type") {
        SECTION("should be a string") {
            REQUIRE(base::determine_string_data_type("") == DataType::STRING);
            REQUIRE(base::determine_string_data_type(" ") == DataType::STRING);
            REQUIRE(base::determine_string_data_type("a") == DataType::STRING);
            // two dots in a row
            REQUIRE(base::determine_string_data_type("1..0") == DataType::STRING);
            //char at the end
            REQUIRE(base::determine_string_data_type("1.0x") == DataType::STRING);
            // plus sign in the middle
            REQUIRE(base::determine_string_data_type("2+5") == DataType::STRING);

        }
        // should be ints
        REQUIRE(base::determine_string_data_type("1") == DataType::INT);

        // should be floats
        REQUIRE(base::determine_string_data_type("1.0") == DataType::FLOAT);


    }

}