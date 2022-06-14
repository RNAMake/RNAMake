//
//  basic_io.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 3/1/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "util/basic_io.hpp"


namespace util {

String
points_to_pdb_str(
        math::Vector3s const & points) {

    String s;
    int i = 1;
    for (auto const & p : points) {
        char buffer[200];
        std::sprintf(buffer, "ATOM %6d  P   C   A   1 %11.3f%8.3f%8.3f  1.00 62.18           P\n", i,
                     p.get_x(), p.get_y(), p.get_z());
        s += String(buffer);
        i++;
    }

    return s;
}


void
points_to_pdb(
        String const & filename,
        math::Vector3s const & points) {

    auto s = points_to_pdb_str(points);
    std::ofstream out;
    out.open(filename);
    out << s;
    out.close();
}

}