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
    SECTION("test determine_string_contents") {
        SECTION("should be a string") {
            REQUIRE(base::determine_string_contents("") == base::StringContents::STRING);
            REQUIRE(base::determine_string_contents(" ") == base::StringContents::STRING);
            REQUIRE(base::determine_string_contents("a") == base::StringContents::STRING);
            // two dots in a row
            REQUIRE(base::determine_string_contents("1..0") == base::StringContents::STRING);
            //char at the end
            REQUIRE(base::determine_string_contents("1.0x") == base::StringContents::STRING);
            // plus sign in the middle
            REQUIRE(base::determine_string_contents("2+5") == base::StringContents::STRING);

        }
        // should be ints
        REQUIRE(base::determine_string_contents("1") == base::StringContents::INT);

        // should be floats
        REQUIRE(base::determine_string_contents("1.0") == base::StringContents::FLOAT);


    }

}