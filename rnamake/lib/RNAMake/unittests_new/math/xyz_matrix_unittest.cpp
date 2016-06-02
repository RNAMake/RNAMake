
#include "../common.hpp"

#include "math/xyz_matrix.h"


TEST_CASE( "Test Matrix math ", "[XYZMatrix]" ) {

    SECTION("Single known test of unitarize compared to python") {
    
        String line;
        std::ifstream input;
        //input.open();
        while ( input.good() ) {
            getline(input, line);
            break;
        }
        input.close();

        
        auto m = Matrix(1.0, 2.0, 3.0,
                        4.0, 5.0, 6.0,
                        1.0, 1.0, 1.0);
    
        auto unit = m.get_unitarize();
    
        std::cout << unit << std::endl;
    
    }
}