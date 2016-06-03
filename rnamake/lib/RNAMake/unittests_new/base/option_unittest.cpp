
#include "../common.hpp"

#include "base/option.h"


TEST_CASE( "Test Options for storing options for classes", "[Options]" ) {
    
    SECTION("Test creation of Option, make sure type is inforced") {
        
        auto opt = Option("test", 6, OptionType::FLOAT);
        REQUIRE(opt.get_float() == 6);
        
        opt = Option("test", "test", OptionType::STRING);
        REQUIRE(opt.get_string() == "test");
        
        opt = Option("test", false, OptionType::BOOL);
        REQUIRE(opt.get_bool() == false);
        
        SECTION("specifying wrong option type should return an error") {
            REQUIRE_THROWS_AS(Option("test", 6.0f, OptionType::STRING), OptionException);
            REQUIRE_THROWS_AS(Option("test", 6, OptionType::STRING), OptionException);
            REQUIRE_THROWS_AS(Option("test", "test", OptionType::FLOAT), OptionException);
            REQUIRE_THROWS_AS(Option("test", "test", OptionType::BOOL), OptionException);
            REQUIRE_THROWS_AS(Option("test", false, OptionType::STRING), OptionException);

        }
        
    }
    
    SECTION("Test getting correct value from option ") {
        
        SECTION("float option cannot be converted to string or bool") {
            auto opt = Option("test", 6, OptionType::FLOAT);
        
            REQUIRE_THROWS_AS(opt.get_string(), OptionException);
            REQUIRE_THROWS_AS(opt.get_bool(), OptionException);
        }

        SECTION("string option cannot be converted to anything else") {
        
            auto opt = Option("test", "test", OptionType::STRING);
            REQUIRE_THROWS_AS(opt.get_int(), OptionException);
            REQUIRE_THROWS_AS(opt.get_float(), OptionException);
            REQUIRE_THROWS_AS(opt.get_bool(), OptionException);

        }
        
        SECTION("int option cannot be converted to string or bool") {
            auto opt = Option("test", 6, OptionType::INT);
            
            REQUIRE_THROWS_AS(opt.get_string(), OptionException);
            REQUIRE_THROWS_AS(opt.get_bool(), OptionException);

        }
        
        
        SECTION("bool option cannot be converted to anything else") {
            
            auto opt = Option("test", false, OptionType::BOOL);
            REQUIRE_THROWS_AS(opt.get_int(), OptionException);
            REQUIRE_THROWS_AS(opt.get_float(), OptionException);
            REQUIRE_THROWS_AS(opt.get_string(), OptionException);
            
        }
    }
    
    SECTION("Test setting new option values") {
        
        SECTION("test setting float options should be able to accept float and int") {
        
            auto opt = Option("test", 6, OptionType::FLOAT);
            opt.value(3);
            REQUIRE(opt.get_float() == 3);
            
            opt.value(3.5f);
            REQUIRE(opt.get_float() == 3.5);
            
            REQUIRE_THROWS_AS(opt.value("test"), OptionException);
            REQUIRE_THROWS_AS(opt.value(false), OptionException);

            
        }
        
        SECTION("test setting int options should be able to accept float and int") {
            auto opt = Option("test", 6, OptionType::INT);
            opt.value(3);
            REQUIRE(opt.get_float() == 3);
            
            opt.value(3.5f);
            REQUIRE(opt.get_float() == 3.0);
            
            REQUIRE_THROWS_AS(opt.value("test"), OptionException);
            REQUIRE_THROWS_AS(opt.value(false), OptionException);

        }
        
        SECTION("test setting string options should be able to accept only str") {
            auto opt = Option("test", "test", OptionType::STRING);
            opt.value("test2");
            REQUIRE(opt.get_string() == "test2");
            
            REQUIRE_THROWS_AS(opt.value(5), OptionException);
            REQUIRE_THROWS_AS(opt.value(5.0f), OptionException);
            REQUIRE_THROWS_AS(opt.value(false), OptionException);
    
        }
        
        SECTION("test setting bool options should be able to accept only bool") {
            auto opt = Option("test", false, OptionType::BOOL);
            opt.value(true);
            REQUIRE(opt.get_bool() == true);
            
            REQUIRE_THROWS_AS(opt.value(5), OptionException);
            REQUIRE_THROWS_AS(opt.value(5.0f), OptionException);
            REQUIRE_THROWS_AS(opt.value("test"), OptionException);
            
        }
        
    }
    
    SECTION("Test getting option objects in options class") {
        auto opts = Options();
        
        opts.add_option("test", "test", OptionType::STRING);
        opts.add_option("test2", 5, OptionType::FLOAT);
        opts.add_option("test3", false, OptionType::BOOL);

        REQUIRE(opts.get_bool("test3") == false);
        
        opts.set_value("test3", true);
        
        REQUIRE(opts.get_bool("test3") == true);
        REQUIRE(opts.has_option("test") == true);
        
        auto names = String();
        for(auto const & opt : opts) {
            names += opt->name() + " ";
        }
        
        std::cout << names << std::endl;
        auto expected = String("test test2 test3 ");
        REQUIRE(names == expected);

    }

}