//
//  main.cpp
//  vienna_clone
//
//  Created by Joseph Yesselman on 4/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "types.h"
#include "vienna.h"
#include "FileIO.h"

Strings
get_lines_from_file(String const fname) {
    String line;
    Strings lines;
    std::ifstream input;
    input.open(fname);
    while ( input.good() ) {
        getline(input, line);
        if( line.length() < 2 ) { break; }
        lines.push_back(line);
        
    }
    return lines;
}

int test_folding() {
    Strings lines = get_lines_from_file("test_folding.dat");
    for(auto const & l : lines) {
        Strings spl = split_str_by_delimiter(l, " ");
        Vienna vc;
        vc.init_fold(1000);
        float energy = vc.fold(spl[0]);
        String structure = vc.get_structure();
        //std::cout << energy << " " << spl[1] << std::endl;
        //std::cout << structure << " " << spl[2] << std::endl;
    }
    return 0;
}

int test_folding_no_reset() {
    Strings lines = get_lines_from_file("test_folding.dat");
    Vienna vc;
    vc.init_fold(1000);
    //for(int i = 0; i < 10000; i++) {
    for(auto const & l : lines) {
        Strings spl = split_str_by_delimiter(l, " ");
        
        float energy = vc.fold(spl[0]);
        String structure = vc.get_structure();
        //std::cout << energy << " " << spl[1] << std::endl;
        //std::cout << structure << " " << spl[2] << std::endl;
    }
    //}

    return 0;
}



int main(int argc, const char * argv[]) {
    //test_folding();
    //test_folding_no_reset();
    
    Vienna v;
    v.dotplot("AAAAAAAAAAUUUUUUUUUU");
    
    return 0;
}
