//
// Created by Hassan Abdelsamad on 4/18/21.
//

#include <CLI/CLI.hpp>
#include "base/backtrace.h"
#include "base/log.h"
#include "design_rna_scaffold/design_rna_scaffold.h"
#include "design_rna_scaffold/design_rna_scaffold.cc"
#include <data_structure/graph_base.h>
#include "tools.h"
#include <base/log.h>


#include "common.hpp"


TEST_CASE ("Test TTR") {

    auto args = get_arguments("TTR");
    auto expected_path = get_expected_path("TTR");
    mock_main(args);

    SUBCASE("Compare logs") {
        CHECK(check_logs(expected_path, "logs.csv"));
    }
}

