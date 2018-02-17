//
// Created by Joseph Yesselman on 2/16/18.
//


#include "../common.hpp"

#include "base/settings.h"
#include "base/file_io.h"
#include "math/xyz_matrix.h"
#include "math/hashing.h"

TEST_CASE( "Test Matrix math ", "[XYZMatrix]" ) {

    auto lower = Point( 12.5, 16.25, 4.25 );
    auto upper = Point( 15.5, 20, 8.5 );
    auto bb = BoundingBox(lower, upper);

    Real6 binwidths;
    binwidths[ 0 ] = binwidths[ 1 ] = binwidths[ 2 ] = 0.25;
    binwidths[ 3 ] = binwidths[ 4 ] = binwidths[ 5 ] = 10;

    Size3 euler_offsets = {0, 0, 0};

    auto binner = SixDCoordinateBinner(bb, binwidths, euler_offsets);

    Real6 pA;
    pA[0] = 13.6; pA[1] = 19.4; pA[2] = 5.3;
    pA[3] = 50;   pA[4] = 127;  pA[5] = 76;

    auto bin = binner.bin6(pA);
    auto bin_center = binner.bin_center_point(bin);

    /*for(auto const & v : bin_center) {
        std::cout << v << std::endl;
    }*/



}