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


 #include "common.hpp"

 TEST_CASE ("Test TTR_supplied_motif") {
//     auto base_path = base::base_dir() + "/apps/simulate_tectos/resources/";
//     resources::Manager::instance().add_motif(base_path + "GAAA_tetraloop");

     // Get the args
     auto base = std::getenv("RNAMAKE");
     std::string base_str(base);
     std::string args_path = base_str + "/integration_tests/TTR_SUPPLIED_MOTIF/ARGS.txt";
     std::ifstream args_file(args_path);
     if (!args_file){
                 FAIL("File does not exist, Make sure the ARGS.txt file exists in the test folder");
     }
     std::string args((std::istreambuf_iterator<char>(args_file)),
                      std::istreambuf_iterator<char>());


     //Get expected logs
     std::string expected_path = base_str + "/integration_tests/TTR_SUPPLIED_MOTIF/EXPECTED.csv";
     std::ifstream expected_file(expected_path);
     if (!expected_file){
                 FAIL("File does not exist, Make sure the EXPECTED.csv file exists in the test folder");
     }

     std::set_terminate(base::print_backtrace);
     auto app = DesignRNAScaffold();

     app.setup_options();
     base::init_logging(app.log_level());
     app.app_.parse(args);
     app.run();

     SUBCASE ("Compare logs") {
         CHECK(check_logs(expected_path, "Logs.txt"));
     }
}


