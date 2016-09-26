
#include "../common.hpp"
#include "../command_line_args.hpp"

#include "base/cl_option.h"
#include "base/string.h"
#include "base/application.hpp"


class TestSearch {
public:
    TestSearch() {
        options_ = Options();
        options_.add_option("accept_score", 10, OptionType::FLOAT);
    }
    
    ~TestSearch() {}
    
public:
    inline
    Options &
    options() { return options_; }
    
private:
    Options options_;
};

class TestApplication : public Application {
public:
    TestApplication() : Application(),
    search_(TestSearch())
    { setup_options(); }
    
    ~TestApplication() {}
    
public:
    
    void
    setup_options() {
        add_option("test", String(""), OptionType::STRING, false);
        add_cl_options(search_.options(), "search");
    }
    
    void
    parse_command_line(
        int argc,
        const char ** argv) {
        
        Application::parse_command_line(argc, argv);
        cl_parser_.assign_options(cl_options_, search_.options(), "search");
        
    }
    
    void
    run() { }
    
public:
    TestSearch search_;
    
};


TEST_CASE( "Test Application, a class to inherit for excutables", "[Application]" ) {

    SECTION("Test parsing options") {
        
        auto cla = CommandLineArgs("-test test2 -search.accept_score 5");
        
        auto app = TestApplication();
        app.parse_command_line(cla.argc, cla.argv());
        
        REQUIRE(app.get_string_option("test") == "test2");
        
        SECTION("test parsing options into classes inside application") {
            REQUIRE(app.search_.options().get_float("accept_score") == 5);
        }
        
    }
    
    

}










