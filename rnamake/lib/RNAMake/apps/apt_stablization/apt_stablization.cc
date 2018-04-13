//
// Created by Joseph Yesselman on 4/12/18.
//

#include "base/backtrace.hpp"
#include "resources/resource_manager.h"
#include "apt_stablization/apt_stablization.h"



APTStablization::APTStablization() : Application() {}

void
APTStablization::setup_options() {
    add_option("pdb", String(""), OptionType::STRING, true);

    // general options
    add_option("out_file", "default.out", OptionType::STRING, false);
    add_option("score_file", "default.scores", OptionType::STRING, false);
    add_option("designs", 1, OptionType::INT, false);

}

void
APTStablization::parse_command_line(
        int argc,
        const char **argv) {

    Application::parse_command_line(argc, argv);
}

void
APTStablization::run() {

    // add motif to resource manager
    RM::instance().add_motif(get_string_option("pdb"), "aptamer", MotifType::TWOWAY);
    std::cout << "APT STABLIZATION: loaded pdb from file: " << get_string_option("pdb") << std::endl;

    auto m = RM::instance().motif("aptamer");
    for(auto const & end : m->ends()) {
        std::cout << end->name() << std::endl;
    }

}


int main(int argc, const char *argv[]) {
    //must add this for all apps!
    std::set_terminate(print_backtrace);

    auto app = APTStablization();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}
