
#include "../common.hpp"
#include "../command_line_args.hpp"

#include "base/cl_option.h"
#include "base/string.h"


TEST_CASE( "Test collecting command line options", "[CLOptions]" ) {
    
    SECTION("test adding new cl options to be parsed") {
        
        auto cl_opts = base::CommandLineOptions();
        cl_opts.add_option("test", "vtest", base::OptionType::STRING);
        
        REQUIRE(cl_opts.get_string("test") == "vtest");
        
    }
    
    //this isnt working out in the real world need to rethink it
    SECTION("boolean commandline options must start at false") {
        //REQUIRE_THROWS_AS(CommandLineOption("test", true, base::OptionType::BOOL, false),
        //                  base::CommandLineOptionException);
    }
    
    SECTION("test parsing command line") {
        
        SECTION("test single option parse") {
            auto cl_opts = base::CommandLineOptions();
            cl_opts.add_option("test", 5, base::OptionType::FLOAT);
            cl_opts.add_option("test_2", "test_3", base::OptionType::STRING);
        
            auto cla = CommandLineArgs("-test 7.0");
            cl_opts.parse_command_line(cla.argc, cla.argv());
        
            REQUIRE(cl_opts.get_float("test") == 7.0);
            REQUIRE(cl_opts.get_string("test_2") == "test_3");
            REQUIRE(cl_opts.is_filled("test") == true);
        }
        
        SECTION("test multi option parse") {
            auto cl_opts = base::CommandLineOptions();
            cl_opts.add_option("test", 5, base::OptionType::FLOAT);
            cl_opts.add_option("test_2", "test_3", base::OptionType::STRING);
            
            auto cla = CommandLineArgs("-test 7.0 -test_2 found");
            cl_opts.parse_command_line(cla.argc, cla.argv());

            REQUIRE(cl_opts.get_float("test") == 7.0);
            REQUIRE(cl_opts.get_string("test_2") == "found");

            
        }
        
        SECTION("test boolean option parse") {
            auto cl_opts = base::CommandLineOptions();
            cl_opts.add_option("test", 5, base::OptionType::FLOAT);
            cl_opts.add_option("test_2", false, base::OptionType::BOOL);
            
            auto cla = CommandLineArgs("-test 7.0 -test_2");
            cl_opts.parse_command_line(cla.argc, cla.argv());

            REQUIRE(cl_opts.get_float("test") == 7);
            REQUIRE(cl_opts.get_bool("test_2") == true);
            

        }
        
        

    }
    
    SECTION("test supplying incorrect value type at command line") {
        
        SECTION("only floats should be accepted for float options") {
        
            auto cl_opts = base::CommandLineOptions();
            cl_opts.add_option("test", 5, base::OptionType::FLOAT);
            cl_opts.add_option("test_2", "test_3", base::OptionType::STRING);
        
            auto cla = CommandLineArgs("-test test");
            REQUIRE_THROWS_AS(cl_opts.parse_command_line(cla.argc, cla.argv()),
                              base::CommandLineOptionException);
        
            cla = CommandLineArgs("-test");
            REQUIRE_THROWS_AS(cl_opts.parse_command_line(cla.argc, cla.argv()),
                              base::CommandLineOptionException);
        
        }
        
        SECTION("only ints should be accepted for int options") {
            auto cl_opts = base::CommandLineOptions();
            cl_opts.add_option("test", 5, base::OptionType::INT);
            
            auto cla = CommandLineArgs("-test test");
            REQUIRE_THROWS_AS(cl_opts.parse_command_line(cla.argc, cla.argv()),
                              base::CommandLineOptionException);
        }
        
        //phasing these tests out as they sc
        SECTION("only string should be accepted for string options") {
            auto cl_opts = base::CommandLineOptions();
            cl_opts.add_option("test", "test", base::OptionType::STRING);
            
            auto cla = CommandLineArgs("-test 5");
            //REQUIRE_THROWS_AS(cl_opts.parse_command_line(cla.argc, cla.argv()),
            //                  base::CommandLineOptionException);
        
            cla = CommandLineArgs("-test 5.0");
            //REQUIRE_THROWS_AS(cl_opts.parse_command_line(cla.argc, cla.argv()),
            //                  base::CommandLineOptionException);
        }
        
        
        SECTION("only bool should be accepted for bool options") {
            auto cl_opts = base::CommandLineOptions();
            cl_opts.add_option("test", false, base::OptionType::BOOL);
            
            auto cla = CommandLineArgs("-test test");
            REQUIRE_THROWS_AS(cl_opts.parse_command_line(cla.argc, cla.argv()),
                              base::CommandLineOptionException);
            
            cla = CommandLineArgs("-test 5");
            REQUIRE_THROWS_AS(cl_opts.parse_command_line(cla.argc, cla.argv()),
                              base::CommandLineOptionException);
            
        }
        
        
    }
    
    SECTION("catch supplying the same argument twice") {
        auto cl_opts = base::CommandLineOptions();
        cl_opts.add_option("test", 5, base::OptionType::FLOAT);
        cl_opts.add_option("test_2", "test_3", base::OptionType::STRING);
        
        auto cla = CommandLineArgs("-test 4 -test_2 -test 10");
        
        REQUIRE_THROWS_AS(cl_opts.parse_command_line(cla.argc, cla.argv()),
                          base::CommandLineOptionException);

    }
    
    SECTION("test add options from option class not by add_option") {
        auto opts = base::Options();
        opts.add_option("test_2", "testv", base::OptionType::STRING);
        opts.add_option("test", 5, base::OptionType::FLOAT);
        
        auto cl_opts = base::CommandLineOptions();
        cl_opts.add_options(opts);
        
        auto cla = CommandLineArgs("-test 7.0 -test_2 found");
        cl_opts.parse_command_line(cla.argc, cla.argv());
        
        REQUIRE(cl_opts.get_float("test") == 7.0);
        REQUIRE(cl_opts.get_string("test_2") == "found");
    }
    
    SECTION("catch name collisions cant add an option that already exists") {
        auto cl_opts = base::CommandLineOptions();
        cl_opts.add_option("test", 5, base::OptionType::FLOAT);
    
        REQUIRE_THROWS_AS(cl_opts.add_option("test", 5, base::OptionType::FLOAT),
                          base::CommandLineOptionException);
    }
    
}



































