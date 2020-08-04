
#include "../common.hpp"

#include "base/option.h"


TEST_CASE( "Test Options for storing options for classes", "[Options]" ) {
    
    SECTION("Test creation of Option, make sure type is inforced") {
        
        auto opt = base::Option("test", 6, base::OptionType::FLOAT);
        REQUIRE(opt.get_float() == 6);
        
        opt = base::Option("test", "test", base::OptionType::STRING);
        REQUIRE(opt.get_string() == "test");
        
        opt = base::Option("test", false, base::OptionType::BOOL);
        REQUIRE(opt.get_bool() == false);
        
        SECTION("specifying wrong option type should return an error") {
            REQUIRE_THROWS_AS(base::Option("test", 6.0f, base::OptionType::STRING),  base::OptionException);
            REQUIRE_THROWS_AS(base::Option("test", 6, base::OptionType::STRING),     base::OptionException);
            REQUIRE_THROWS_AS(base::Option("test", "test", base::OptionType::FLOAT), base::OptionException);
            REQUIRE_THROWS_AS(base::Option("test", "test", base::OptionType::BOOL),  base::OptionException);
            REQUIRE_THROWS_AS(base::Option("test", false, base::OptionType::STRING), base::OptionException);

        }
        
    }
    
    SECTION("Test getting correct value from option ") {
        
        SECTION("float option cannot be converted to string or bool") {
            auto opt = base::Option("test", 6, base::OptionType::FLOAT);
        
            REQUIRE_THROWS_AS(opt.get_string(), base::OptionException);
            REQUIRE_THROWS_AS(opt.get_bool(), base::OptionException);
        }

        SECTION("string option cannot be converted to anything else") {
        
            auto opt = base::Option("test", "test", base::OptionType::STRING);
            REQUIRE_THROWS_AS(opt.get_int(),   base::OptionException);
            REQUIRE_THROWS_AS(opt.get_float(), base::OptionException);
            REQUIRE_THROWS_AS(opt.get_bool(),  base::OptionException);

        }
        
        SECTION("int option cannot be converted to string or bool") {
            auto opt = base::Option("test", 6, base::OptionType::INT);
            
            REQUIRE_THROWS_AS(opt.get_string(), base::OptionException);
            REQUIRE_THROWS_AS(opt.get_bool(),   base::OptionException);

        }
        
        
        SECTION("bool option cannot be converted to anything else") {
            
            auto opt = base::Option("test", false, base::OptionType::BOOL);
            REQUIRE_THROWS_AS(opt.get_int(),    base::OptionException);
            REQUIRE_THROWS_AS(opt.get_float(),  base::OptionException);
            REQUIRE_THROWS_AS(opt.get_string(), base::OptionException);
            
        }
    }
    
    SECTION("Test setting new option values") {
        
        SECTION("test setting float options should be able to accept float and int") {
        
            auto opt = base::Option("test", 6, base::OptionType::FLOAT);
            opt.value(3);
            REQUIRE(opt.get_float() == 3);
            
            opt.value(3.5f);
            REQUIRE(opt.get_float() == 3.5);
            
            REQUIRE_THROWS_AS(opt.value("test"), base::OptionException);
            REQUIRE_THROWS_AS(opt.value(false),  base::OptionException);

            
        }
        
        SECTION("test setting int options should be able to accept float and int") {
            auto opt = base::Option("test", 6, base::OptionType::INT);
            opt.value(3);
            REQUIRE(opt.get_float() == 3);
            
            opt.value(3.5f);
            REQUIRE(opt.get_float() == 3.0);
            
            REQUIRE_THROWS_AS(opt.value("test"), base::OptionException);
            REQUIRE_THROWS_AS(opt.value(false),  base::OptionException);

        }
        
        SECTION("test setting string options should be able to accept only str") {
            auto opt = base::Option("test", "test", base::OptionType::STRING);
            opt.value("test2");
            REQUIRE(opt.get_string() == "test2");

            REQUIRE_THROWS_AS(opt.value(5),     base::OptionException);
            REQUIRE_THROWS_AS(opt.value(5.0f),  base::OptionException);
            REQUIRE_THROWS_AS(opt.value(false), base::OptionException);
    
        }
        
        SECTION("test setting bool options should be able to accept only bool") {
            auto opt = base::Option("test", false, base::OptionType::BOOL);
            opt.value(true);
            REQUIRE(opt.get_bool() == true);
            
            REQUIRE_THROWS_AS(opt.value(5),      base::OptionException);
            REQUIRE_THROWS_AS(opt.value(5.0f),   base::OptionException);
            REQUIRE_THROWS_AS(opt.value("test"), base::OptionException);
            
        }
        
    }
    
    SECTION("Test getting option objects in options class") {
        auto opts = base::Options();
        
        opts.add_option("test", "test", base::OptionType::STRING);
        opts.add_option("test2", 5, base::OptionType::FLOAT);
        opts.add_option("test3", false, base::OptionType::BOOL);

        REQUIRE(opts.get_bool("test3") == false);
        
        opts.set_value("test3", true);
        
        REQUIRE(opts.get_bool("test3") == true);
        REQUIRE(opts.has_option("test") == true);
        
        auto names = String();
        for(auto const & opt : opts) {
            names += opt->name() + " ";
        }
        
        auto expected = String("test test2 test3 ");
        REQUIRE(names == expected);

    }

}