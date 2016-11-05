
#include "../common.hpp"
#include "../command_line_args.hpp"


TEST_CASE( "Test DesignRNA Application", "[DesignRNA]" ) {
    
    SECTION("test using pdbs") {
        auto cla = CommandLineArgs("-pdb test.pdb");
    }
    
    
}


