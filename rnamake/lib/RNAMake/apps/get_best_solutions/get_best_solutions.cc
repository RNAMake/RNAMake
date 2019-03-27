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
    add_option("info_file", String(""), base::OptionType::STRING, false);
    add_option("new_out_file", String(""), base::OptionType::STRING, false);
    add_option("sequence_file", String(""), base::OptionType::STRING, false);

}

void
GetBestSolutions::parse_command_line(
        int argc,
        const char ** argv) {
    base::Application::parse_command_line(argc, argv);

    parameters_.out_file = get_string_option("out_file");
    parameters_.info_file = get_string_option("info_file");
    parameters_.new_out_file = get_string_option("new_out_file");
    parameters_.sequence_file = get_string_option("sequence_file");
    parameters_.rows = get_string_option("rows");

    // parse row selection info
    if(parameters_.rows != "") {
        auto spl = base::split_str_by_delimiter(parameters_.rows, ",");
        for(auto const & s : spl) {
            auto row = std::stoi(s);
            rows_[row] = 1;
        }
        parameters_.using_rows = true;
    }

    // parse sequences
    if(parameters_.sequence_file != "") {
        auto lines = base::get_lines_from_file(get_string_option("sequence_file"));
        for(auto const & l : lines) {
            sequences_.push_back(l);
        }
        parameters_.using_sequence_file = true;
    }

    // are we printing out new files
    if(parameters_.info_file != "") {
        info_file_ = std::make_shared<std::ofstream>();
        info_file_->open(parameters_.info_file);
        *info_file_ << "design_num,sequence,structure" << std::endl;
        parameters_.using_info_file = true;
    }

    if(parameters_.new_out_file != "") {
        new_out_file_ = std::make_shared<std::ofstream>();
        new_out_file_->open(parameters_.new_out_file);
        parameters_.using_new_out_file = true;
    }


}

void
GetBestSolutions::run() {

    auto lines = base::get_lines_from_file(get_string_option("out_file"));
    auto i = -1;
    for(auto const & l : lines) {
        i += 1;
        // row was not selected
        if(parameters_.using_rows && rows_.find(i) == rows_.end()) { continue; }
        auto mg = std::make_shared<motif_data_structure::MotifGraph>(l, motif_data_structure::MotifGraphStringType::MG);
        _fix_flex_helix_type(*mg);
        _process_motif_graph(*mg, i);
    }

    if(parameters_.using_info_file) { info_file_->close(); }
    if(parameters_.using_new_out_file) { new_out_file_->close(); }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void
GetBestSolutions::_process_motif_graph(
        motif_data_structure::MotifGraph & mg,
        int row) {

    if(parameters_.using_sequence_file) {
        mg.replace_ideal_helices();
        mg.replace_helical_sequence(sequences_[row]);
    }

    if(parameters_.using_new_out_file) {
        *new_out_file_ << mg.to_str() << std::endl;
    }

    if(parameters_.using_info_file) {
        *info_file_ << row << "," << mg.sequence() << " " << mg.dot_bracket() << std::endl;
    }

    mg.to_pdb("test." + std::to_string(row) + ".pdb", 1, 1);

}

void
GetBestSolutions::_fix_flex_helix_type(
        motif_data_structure::MotifGraph & mg) {
    for(auto & n : mg) {
        if(n->data()->name().substr(0, 5) == "HELIX") {
            n->data()->mtype(util::MotifType::HELIX);
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