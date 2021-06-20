 //
 // Created by Hassan Abdelsamad on 4/18/21.
 //

 #include <CLI/CLI.hpp>

 #include "design_rna_scaffold/design_rna_scaffold.cpp"
 #include "tools.h"


 #include "common.hpp"

 TEST_CASE ("Test TTR_supplied_motif") {

     auto args = get_arguments("TTR_SUPPLIED_MOTIF");
     auto expected_path = get_expected_path("TTR_SUPPLIED_MOTIF");
     mock_main(args);

     SUBCASE ("Compare logs") {
         CHECK(check_logs(expected_path, "logs.csv"));
     }
}


