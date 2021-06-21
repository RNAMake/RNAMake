//
// Created by Hassan Abdelsamad on 4/30/21.
//


#include "design_rna_scaffold/design_rna_scaffold.cpp"
#include "tools.h"
#include "common.hpp"


TEST_CASE ("Test rebuild_p4_p6") {

    auto args = get_arguments("REBUILD_P4_P6");
    auto expected_path = get_expected_path("REBUILD_P4_P6");
    mock_main(args);

    SUBCASE("Compare logs") {
        CHECK(check_logs(expected_path, "logs.csv"));
    }
}