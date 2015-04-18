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

    return 0;
}

int test_bp_probs() {
    Strings lines = get_lines_from_file("test_bp_probs.dat");
    plists bp_probs;
    int i, j;
    float prob;
    int count = 0;
    Vienna v;
    for(auto const & l : lines) {
        Strings spl = split_str_by_delimiter(l, "|");
        String seq = spl[0];
        std::cout << seq << std::endl;
        bp_probs = v.bp_probabilities(seq);
        //std::cout << count << " " << seq << std::endl;
        std::cout << v.get_structure() << std::endl;
        for(int k = 1; k < spl.size(); k++) {
            Strings e = split_str_by_delimiter(spl[k], " ");
            i = stoi(e[0]);
            j = stoi(e[1]);
            prob = stof(e[2]);
            for(auto const & p : bp_probs) {
                if(p.i == 0 && p.j == 0) { break; }
                if(i == p.i && j == p.j) {
                    if(fabs(prob - p.p) > 0.1) {
                        std::cout << count << " " << seq << " " << i << " " << j << " " << prob << " " << p.p << std::endl;
                    }
                    break;
                }
            }
            
        }
        count++;
    
    }
    

    return 0;
}

int test_memory_leak() {
    String ful_seq = "GGGGAUAUGGGGGGGGGGGGGGAUGGAAGGGGGGGGGGGGGGGGGGGGGCAACAGCGAGGGAGAGGGAAACCAAGUUCCCAUCGACAUGCCCCCCCCCCCCCCCCCCCCUGGACCCCCCCCCCCCCCUAAGUCCCCGGGGAUAUGGGGGGGGGGGGGGAUGGAAGGGGGGGGGGGGGGGGGGGGGCAACAGCGAGGGAGAGGGAAACCAAGUUCCCAUCGACAUGCCCCCCCCCCCCCCCCCCCCUGGACCCCCCCCCCCCCCUAAGUCCCC";
    int count = 0;
    String seq;
    Vienna v;

    for(int i = 0; i < 10000; i++) {
        seq = ful_seq.substr(0,10+count );
        if(count > 240) { count = 0;}
        std::cout << i << " " << count << std::endl;
        v.fold(seq);
        count++;
    }

    return 0;
    
}


int main(int argc, const char * argv[]) {
    //test_folding();
    //test_folding_no_reset();
    //test_bp_probs();
    test_memory_leak();
    exit(0);

    return 0;
}
