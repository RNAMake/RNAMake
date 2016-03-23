//
//  basic_io.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 3/1/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "util/basic_io.hpp"


String
points_to_pdb_str(
    Points const & points) {
    
    String s;
    int i = 1;
    for(auto const & p : points) {
        char buffer [200];
        std::sprintf(buffer, "ATOM %6d  P   C   A   1 %11.3f%8.3f%8.3f  1.00 62.18           P\n", i, p.x(), p.y(), p.z());
        s += String(buffer);
        i++;
    }
    
    return s;
}


void
points_to_pdb(
    String const & filename,
    Points const & points) {
    
    auto s = points_to_pdb_str(points);
    std::ofstream out;
    out.open(filename);
    out << s;
    out.close();
}