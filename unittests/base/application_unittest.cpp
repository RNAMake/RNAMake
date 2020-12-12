

#include "../command_line_args.hpp"
#include "../common.hpp"

#include "base/cl_option.h"
#include "base/string.h"
#include "base/application.hpp"


class TestSearch {
public:
    TestSearch() {
        options_ = base::Options();
        options_.add_option("accept_score", 10, base::OptionType::FLOAT);
    }
    
    ~TestSearch() {}
    
public:
    inline
    base::Options &
    options() { return options_; }
    
private:
    base::Options options_;
};

class TestApplication : public base::Application {
public:
    TestApplication() : base::Application(),
    search_(TestSearch())
    { setup_options(); }
    
    ~TestApplication() {}
    
public:
    
    void
    setup_options() {
        add_option("test", String(""), base::OptionType::STRING, false);
        add_cl_options(search_.options(), "search");
    }
    
    void
    parse_command_line(
        int argc,
        const char ** argv) {
        
        base::Application::parse_command_line(argc, argv);
        cl_parser_.assign_options(cl_options_, search_.options(), "search");
        
    }
    
    void
    run() { }
    
public:
    TestSearch search_;
    
};


TEST_CASE( "Test Application, a class to inherit for excutables") {

    SUBCASE("Test parsing options") {
        
        auto cla = CommandLineArgs("-test test2 -search.accept_score 5");
        
        auto app = TestApplication();
        app.parse_command_line(cla.argc, cla.argv());
        
        CHECK(app.get_string_option("test") == "test2");
        
        SUBCASE("test parsing options into classes inside application") {
            CHECK(app.search_.options().get_float("accept_score") == 5);
        }
        
    }
    
    

}










