//
// Created by Joseph Yesselman on 3/19/19.
//

#ifndef RNAMAKE_NEW_GET_BEST_SOLUTIONS_H
#define RNAMAKE_NEW_GET_BEST_SOLUTIONS_H

#include <base/application.hpp>

class Table {
public:
    Table(
            String const & csv_file) {

    }

private:

};


class GetBestSolutions : public base::Application {
public:
    struct Parameters {
        String out_file, info_file, new_out_file, sequence_file, rows;
        bool using_rows, using_info_file, using_new_out_file, using_sequence_file;
    };


public:
    GetBestSolutions();

    ~GetBestSolutions() {}

public:

    void
    setup_options();

    void
    parse_command_line(
            int,
            const char **);

    void
    run();

private:
    void
    _process_motif_graph(
            motif_data_structure::MotifGraph &,
            int);

    void
    _fix_flex_helix_type(
            motif_data_structure::MotifGraph &);

private:
    Parameters parameters_;
    Strings sequences_;
    std::map<int, int> rows_;
    std::shared_ptr<std::ofstream> info_file_, new_out_file_;

};


#endif //RNAMAKE_NEW_GET_BEST_SOLUTIONS_H
