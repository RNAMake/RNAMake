//
//  eternabot.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "eternabot.h"

#include <base/backtrace.hpp>
#include <base/log.h>
#include "base/cl_option.h"
#include "secondary_structure/secondary_structure_parser.h"
#include "eternabot/sequence_designer.h"

void
EternabotApp::setup_options() {
    add_option("seq", String(""), base::OptionType::STRING, true);
    add_option("ss", String(""), base::OptionType::STRING, true);
    add_option("steps", 100, base::OptionType::INT);
    add_option("n", 1, base::OptionType::INT);

}

void
EternabotApp::parse_command_line(
        int argc,
        const char ** argv) {
    base::Application::parse_command_line(argc, argv);
    parameters_.seq   = get_string_option("seq");
    parameters_.ss    = get_string_option("ss");
    parameters_.steps = get_int_option("steps");
    parameters_.n     = get_int_option("n");
}

void
EternabotApp::run() {
    base::init_logging();

    auto designer = eternabot::SequenceDesigner();
    designer.set_option_value("steps", parameters_.steps);
    designer.setup();

    auto parser = secondary_structure::Parser();
    for(int i = 0; i < parameters_.n; i++) {
        auto p = parser.parse_to_pose(parameters_.seq, parameters_.ss);
        auto results = designer.design(p);
        std::cout << results[0]->score << " " << results[0]->sequence << " " << p->dot_bracket() << std::endl;

    }
}


int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(base::print_backtrace);

    auto app = EternabotApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;


}
