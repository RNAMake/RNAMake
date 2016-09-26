
#include "../common.hpp"
#include "../command_line_args.hpp"


TEST_CASE( "Test DesignRNA Application", "[DesignRNA]" ) {
    
    SECTION("test using pdbs") {
        auto = CommandLineArgs("-pdb test.pdb");
    }
    
    
}


