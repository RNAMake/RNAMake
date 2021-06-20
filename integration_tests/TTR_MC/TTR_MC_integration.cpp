 //
 // Created by Hassan Abdelsamad on 4/30/21.
 //

 #include <CLI/CLI.hpp>

 #include "design_rna_scaffold/design_rna_scaffold.cpp"
 #include "tools.h"


 #include "common.hpp"


 TEST_CASE ("Test TTR_MC") {

     auto args = get_arguments("TTR_MC");
     auto expected_path = get_expected_path("TTR_MC");
     mock_main(args);
         SUBCASE("Compare logs") {
             CHECK(check_logs(expected_path, "logs.csv"));
     }
 }