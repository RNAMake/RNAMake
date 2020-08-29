//
// Created by Joseph Yesselman on 2/16/20.
//

//
// Created by Joseph Yesselman on 3/10/19.
//


#include <iostream>
#include "../common.hpp"
#include <base/string.h>
#include <base/types.h>

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

    SECTION("Testing trim methods") {
       
        SECTION("left trim") {
            const auto target = String{"trimmed"};
            auto untrimmed = String{" trimmed"};
            REQUIRE(base::ltrim(untrimmed) == target);
            untrimmed = "      trimmed"; 
            REQUIRE(base::ltrim(untrimmed) == target);
            untrimmed = "\t\ttrimmed"; 
            REQUIRE(base::ltrim(untrimmed) == target);
            untrimmed = "\n\ttrimmed"; 
            REQUIRE(base::ltrim(untrimmed) == target);
            auto empty_str = String{"            "};
            REQUIRE(base::ltrim(empty_str).empty());
        }

        SECTION("right trim") {
            const auto target = String{"trimmed"};
            auto untrimmed = String{"trimmed "};
            REQUIRE(base::rtrim(untrimmed) == target);
            untrimmed = "trimmed            "; 
            REQUIRE(base::rtrim(untrimmed) == target);
            untrimmed = "trimmed\t\t"; 
            REQUIRE(base::rtrim(untrimmed) == target);
            untrimmed = "trimmed\n\t"; 
            REQUIRE(base::rtrim(untrimmed) == target);
            auto empty_str = String{"            "};
            REQUIRE(base::rtrim(empty_str).empty());
        }

        SECTION("trim") {
            const auto target("trimmed");
            auto untrimmed = String{" trimmed "};
            REQUIRE(base::trim(untrimmed) == target);
            untrimmed = "         trimmed            "; 
            REQUIRE(base::trim(untrimmed) == target);
            untrimmed = "\t\ttrimmed\t\t"; 
            REQUIRE(base::trim(untrimmed) == target);
            untrimmed = "\n\ttrimmed\n\t"; 
            REQUIRE(base::trim(untrimmed) == target);
            auto empty_str = String{"            "};
            REQUIRE(base::trim(empty_str).empty());
        }

    }


    SECTION("testing string splitting methods") {
        
        SECTION("leading spaces") {
            // trailing space(s) 
            const auto raw_line = String{"1 2 3 4 5 6 "}; 
            const auto target = Strings{ "1","2","3","4","5","6"};
            const auto actual = base::tokenize_line(raw_line);

            REQUIRE(target.size() == actual.size());
            auto targ_it = target.cbegin();
            auto actual_it = actual.cbegin();
    
            for( ; targ_it != target.cend(); ++targ_it,++actual_it) {
                REQUIRE(*targ_it == *actual_it);
            }
        }
        
        SECTION("spaces on both sides") { 
            // leading space(s) 
            const auto raw_line = String{"     f s _ $ % "}; 
            const auto target = Strings{ "f","s","_","$", "%"};
            const auto actual = base::tokenize_line(raw_line);

            REQUIRE(target.size() == actual.size());
            auto targ_it = target.cbegin();
            auto actual_it = actual.cbegin();
    
            for( ; targ_it != target.cend(); ++targ_it,++actual_it) {
                REQUIRE(*targ_it == *actual_it);
            }
        }
        
        SECTION("empty string") { 
            // leading space(s) 
            const auto raw_line = String{" "}; 
            const auto target = Strings{};
            const auto actual = base::tokenize_line(raw_line);

            REQUIRE(target.size() == actual.size());
            auto targ_it = target.cbegin();
            auto actual_it = actual.cbegin();
    
            for( ; targ_it != target.cend(); ++targ_it,++actual_it) {
                REQUIRE(*targ_it == *actual_it);
            }
        }
        
        SECTION("escape characters") { 
            const auto raw_line = String{"\0\t\t\t\a"}; 
            const auto target = Strings{};
            const auto actual = base::tokenize_line(raw_line);

            REQUIRE(target.size() == actual.size());
            auto targ_it = target.cbegin();
            auto actual_it = actual.cbegin();
    
            for( ; targ_it != target.cend(); ++targ_it,++actual_it) {
                REQUIRE(*targ_it == *actual_it);
            }
        }

        SECTION("quote behavior") { 
            const auto raw_line = String{"\'this has spaces\' last_token"}; 
            const auto target = Strings{"this has spaces", "last_token"};
            const auto actual = base::tokenize_line(raw_line);

            REQUIRE(target.size() == actual.size());
            auto targ_it = target.cbegin();
            auto actual_it = actual.cbegin();
    
            for( ; targ_it != target.cend(); ++targ_it,++actual_it) {
                REQUIRE(*targ_it == *actual_it);
            }
        }
    }

}
