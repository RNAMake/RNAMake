

#include "../common.hpp"

#include "base/option.h"


TEST_CASE( "Test Options for storing options for classes") {
    
    SUBCASE("Test creation of Option, make sure type is inforced") {
        
        auto opt = base::Option("test", 6, base::OptionType::FLOAT);
        CHECK(opt.get_float() == 6);
        
        opt = base::Option("test", "test", base::OptionType::STRING);
        CHECK(opt.get_string() == "test");
        
        opt = base::Option("test", false, base::OptionType::BOOL);
        CHECK(opt.get_bool() == false);
        
        SUBCASE("specifying wrong option type should return an error") {
            REQUIRE_THROWS_AS(base::Option("test", 6.0f, base::OptionType::STRING),  base::OptionException);
            REQUIRE_THROWS_AS(base::Option("test", 6, base::OptionType::STRING),     base::OptionException);
            REQUIRE_THROWS_AS(base::Option("test", "test", base::OptionType::FLOAT), base::OptionException);
            REQUIRE_THROWS_AS(base::Option("test", "test", base::OptionType::BOOL),  base::OptionException);
            REQUIRE_THROWS_AS(base::Option("test", false, base::OptionType::STRING), base::OptionException);

        }
        
    }
    
    SUBCASE("Test getting correct value from option ") {
        
        SUBCASE("float option cannot be converted to string or bool") {
            auto opt = base::Option("test", 6, base::OptionType::FLOAT);
        
            REQUIRE_THROWS_AS(opt.get_string(), base::OptionException);
            REQUIRE_THROWS_AS(opt.get_bool(), base::OptionException);
        }

        SUBCASE("string option cannot be converted to anything else") {
        
            auto opt = base::Option("test", "test", base::OptionType::STRING);
            REQUIRE_THROWS_AS(opt.get_int(),   base::OptionException);
            REQUIRE_THROWS_AS(opt.get_float(), base::OptionException);
            REQUIRE_THROWS_AS(opt.get_bool(),  base::OptionException);

        }
        
        SUBCASE("int option cannot be converted to string or bool") {
            auto opt = base::Option("test", 6, base::OptionType::INT);
            
            REQUIRE_THROWS_AS(opt.get_string(), base::OptionException);
            REQUIRE_THROWS_AS(opt.get_bool(),   base::OptionException);

        }
        
        
        SUBCASE("bool option cannot be converted to anything else") {
            
            auto opt = base::Option("test", false, base::OptionType::BOOL);
            REQUIRE_THROWS_AS(opt.get_int(),    base::OptionException);
            REQUIRE_THROWS_AS(opt.get_float(),  base::OptionException);
            REQUIRE_THROWS_AS(opt.get_string(), base::OptionException);
            
        }
    }
    
    SUBCASE("Test setting new option values") {
        
        SUBCASE("test setting float options should be able to accept float and int") {
        
            auto opt = base::Option("test", 6, base::OptionType::FLOAT);
            opt.value(3);
            CHECK(opt.get_float() == 3);
            
            opt.value(3.5f);
            CHECK(opt.get_float() == 3.5);
            
            REQUIRE_THROWS_AS(opt.value("test"), base::OptionException);
            REQUIRE_THROWS_AS(opt.value(false),  base::OptionException);

            
        }
        
        SUBCASE("test setting int options should be able to accept float and int") {
            auto opt = base::Option("test", 6, base::OptionType::INT);
            opt.value(3);
            CHECK(opt.get_float() == 3);
            
            opt.value(3.5f);
            CHECK(opt.get_float() == 3.0);
            
            REQUIRE_THROWS_AS(opt.value("test"), base::OptionException);
            REQUIRE_THROWS_AS(opt.value(false),  base::OptionException);

        }
        
        SUBCASE("test setting string options should be able to accept only str") {
            auto opt = base::Option("test", "test", base::OptionType::STRING);
            opt.value("test2");
            CHECK(opt.get_string() == "test2");

            REQUIRE_THROWS_AS(opt.value(5),     base::OptionException);
            REQUIRE_THROWS_AS(opt.value(5.0f),  base::OptionException);
            REQUIRE_THROWS_AS(opt.value(false), base::OptionException);
    
        }
        
        SUBCASE("test setting bool options should be able to accept only bool") {
            auto opt = base::Option("test", false, base::OptionType::BOOL);
            opt.value(true);
            CHECK(opt.get_bool() == true);
            
            REQUIRE_THROWS_AS(opt.value(5),      base::OptionException);
            REQUIRE_THROWS_AS(opt.value(5.0f),   base::OptionException);
            REQUIRE_THROWS_AS(opt.value("test"), base::OptionException);
            
        }
        
    }
    
    SUBCASE("Test getting option objects in options class") {
        auto opts = base::Options();
        
        opts.add_option("test", "test", base::OptionType::STRING);
        opts.add_option("test2", 5, base::OptionType::FLOAT);
        opts.add_option("test3", false, base::OptionType::BOOL);

        CHECK(opts.get_bool("test3") == false);
        
        opts.set_value("test3", true);
        
        CHECK(opts.get_bool("test3") == true);
        CHECK(opts.has_option("test") == true);
        
        auto names = String();
        for(auto const & opt : opts) {
            names += opt->name() + " ";
        }
        
        auto expected = String("test test2 test3 ");
        CHECK(names == expected);

    }

}