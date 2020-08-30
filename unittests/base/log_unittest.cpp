//
// Created by Joseph Yesselman on 3/10/19.
//


#include <iostream>
#include "../common.hpp"
#include "base/log.h"

TEST_CASE( "Test logging functions", "[Logger]" ) {
    base::init_logging(base::LogLevel::INFO);

    //LOGI << "Hello log!"; // short macro
    //LOG_ERROR << "Hello log!"; // long macro
    //LOG(plog::debug) << "Hello log!"; // function-style macro

}