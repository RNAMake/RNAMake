//
//  application_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/1/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "application_unittest.hpp"
#include "cl_option_unittest.h"
#include "base/application.hpp"
#include "motif_state_search/motif_state_search.h"



namespace unittests {

    
class TestApplication : public Application {
public:
    TestApplication() : Application(),
    search_(MotifStateSearch())
    {}
    
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
    MotifStateSearch search_;

};


int
ApplicationUnittest::test_creation() {
    auto app = TestApplication();
    return 0;
}

int
ApplicationUnittest::test_setup_cl_options() {
    auto app = TestApplication();
    auto cla = CommandLineArgs("-test test2 -search::accept_score 5");
    
    app.setup_options();
    app.parse_command_line(cla.argc, cla.argv());
    if(app.get_string_option("test") != "test2") {
        throw UnittestException("did not set test variable correctly");
    }
    
    if(fabs(app.search_.get_float_option("accept_score") - 5) > 0.01) {
        throw UnittestException("did not set search option correctly");
    }
    
    
    return 0;
}
    

int
ApplicationUnittest::run() {
    test_creation();
    test_setup_cl_options();
    return 1;
}
    
    
int
ApplicationUnittest::run_all() {
    return 1;
}


}