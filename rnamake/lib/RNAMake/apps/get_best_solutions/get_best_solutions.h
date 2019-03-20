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
};


#endif //RNAMAKE_NEW_GET_BEST_SOLUTIONS_H
