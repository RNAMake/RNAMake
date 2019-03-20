//
// Created by Joseph Yesselman on 3/19/19.
//

#include <base/backtrace.hpp>
#include <base/log.h>
#include <base/file_io.h>
#include <motif_data_structure/motif_graph.h>
#include <get_best_solutions/get_best_solutions.h>

GetBestSolutions::GetBestSolutions() {

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// app functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
GetBestSolutions::setup_options() {
    add_option("out_file", String(""), base::OptionType::STRING, true);
    //add_option("score_file", String(""), base::OptionType::STRING, true);

    // optional
    add_option("rows", String(""), base::OptionType::STRING, false);

}

void
GetBestSolutions::parse_command_line(
        int argc,
        const char ** argv) {
    base::Application::parse_command_line(argc, argv);
}

void
GetBestSolutions::run() {
    auto lines = base::get_lines_from_file(get_string_option("out_file"));

    if(get_string_option("rows") != "") {
        auto spl = base::split_str_by_delimiter(get_string_option("rows"), ",");
        for(auto const & s : spl) {
            auto row = std::stoi(s);
            auto mg = std::make_shared<motif_data_structure::MotifGraph>(lines[row],
                                                                         motif_data_structure::MotifGraphStringType::MG);
            for(auto & n : *mg) {
                if(n->data()->name().substr(0, 5) == "HELIX") {
                    n->data()->mtype(util::MotifType::HELIX);
                }
            }
            mg->to_pdb("test."+s+".pdb", 1, 1);
        }

    }

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int
main(
        int argc,
        const char ** argv) {

    //must add this for all apps!
    std::set_terminate(base::print_backtrace);

    //turn on logging
    base::init_logging();

    auto app = GetBestSolutions();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}