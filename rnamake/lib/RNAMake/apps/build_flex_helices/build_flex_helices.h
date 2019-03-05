//
// Created by Joseph Yesselman on 3/3/19.
//

#ifndef TEST_BUILD_FLEX_HELICES_H
#define TEST_BUILD_FLEX_HELICES_H

#include <stdio.h>
#include "base/application.hpp"
#include "motif/motif.h"

class BuildFlexHelicesAppException : public std::runtime_error {
public:
    BuildFlexHelicesAppException(
            String const & message):
            std::runtime_error(message)
    {}
};


class BuildFlexHelicesApp : public Application {
public:
    BuildFlexHelicesApp();

    ~BuildFlexHelicesApp() {}

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
    MotifStateOPs
    get_motifs_from_seq_and_ss(
            String const &,
            String const &);

    void
    generate_helices(
            int);

    MotifOP
    get_avg_helix(
            int);

    int
    convert_char_to_res_code(
            char);

    bool
    find_seq_violations(
            Ints const &);

    bool
    find_gc_strech(
            Ints const &);

private:
    std::vector<Ints> disallowed_num_sequences_;
};


#endif //TEST_BUILD_FLEX_HELICES_H
