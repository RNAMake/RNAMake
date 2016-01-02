//
//  vienna_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "base/string.h"
#include "util/settings.h"
#include "util/file_io.h"
#include "vienna/vienna.h"

#include "vienna_unittest.h"


namespace unittests {
namespace vienna {

int
ViennaUnittest::test_creation() {
    auto vc = Vienna();
    auto energy = vc.fold("GGGGUUCGCCCC");
    auto structure = vc.get_structure();

    if(!(5.3f < energy < 5.5f)) {
        throw UnittestException("did not get correct energy");
    }
    
    return 0;
}

int
ViennaUnittest::test_folding() {
    String path = lib_path() + "/unittests/resources/vienna/test_folding.dat";
    Strings lines = get_lines_from_file(path);
    for(auto const & l : lines) {
        Strings spl = split_str_by_delimiter(l, " ");
        if(spl.size() < 2) { break; }
        //std::cout << l << std::endl;
        Vienna vc;
        vc.init_fold(1000);
        float energy = vc.fold(spl[0]);
        String structure = vc.get_structure();
        
        if(structure != spl[2]) {
            throw UnittestException("did not produce correct structure when folding");
        }
        
        float diff = energy - std::stof(spl[1]);
        if(-0.1 > diff || diff > 0.1) {
            throw UnittestException("did not get correct energy when folding");
        }
    }
    return 0;
}
    
int
ViennaUnittest::test_folding_no_reset() {
    String path = lib_path() + "/unittests/resources/vienna/test_folding.dat";
    Strings lines = get_lines_from_file(path);
    Vienna vc;
    vc.init_fold(1000);
    for(auto const & l : lines) {
        Strings spl = split_str_by_delimiter(l, " ");
        if(spl.size() < 2) { break; }
        //std::cout << l << std::endl;
        float energy = vc.fold(spl[0]);
        String structure = vc.get_structure();
        
        if(structure != spl[2]) {
            throw UnittestException("did not produce correct structure when folding");
        }
        
        float diff = energy - std::stof(spl[1]);
        if(-0.1 > diff || diff > 0.1) {
            throw UnittestException("did not get correct energy when folding");
        }
    }
    return 0;
}
    
int
ViennaUnittest::test_bp_probs() {
    String path = lib_path() + "/unittests/resources/vienna/test_bp_probs.dat";
    Strings lines = get_lines_from_file(path);
    plists bp_probs;
    int i, j;
    float prob;
    int count = 0;
    Vienna v;
    for(auto const & l : lines) {
        if(l.length() < 5) { break; }
        Strings spl = split_str_by_delimiter(l, "|");
        String seq = spl[0];
        bp_probs = v.bp_probabilities(seq);
        for(int k = 1; k < spl.size(); k++) {
            Strings e = split_str_by_delimiter(spl[k], " ");
            i = stoi(e[0]);
            j = stoi(e[1]);
            prob = stof(e[2]);
            for(auto const & p : bp_probs) {
                if(p.i == 0 && p.j == 0) { break; }
                if(i == p.i && j == p.j) {
                    if(fabs(prob - p.p) > 0.1) {
                        throw UnittestException("did not get correct bp probability");
                    }
                    break;
                }
            }
            
        }
        count++;
        
    }

    return 0;
}
 
int
ViennaUnittest::test_memory_leak() {
    String ful_seq = "GGGGAUAUGGGGGGGGGGGGGGAUGGAAGGGGGGGGGGGGGGGGGGGGGCAACAGCGAGGGAGAGGGAAACCAAGUUCCCAUCGACAUGCCCCCCCCCCCCCCCCCCCCUGGACCCCCCCCCCCCCCUAAGUCCCCGGGGAUAUGGGGGGGGGGGGGGAUGGAAGGGGGGGGGGGGGGGGGGGGGCAACAGCGAGGGAGAGGGAAACCAAGUUCCCAUCGACAUGCCCCCCCCCCCCCCCCCCCCUGGACCCCCCCCCCCCCCUAAGUCCCC";
    int count = 0;
    String seq;
    Vienna v;
    
    for(int i = 0; i < 10000; i++) {
        seq = ful_seq.substr(0,10+count );
        if(count > 240) { count = 0;}
        v.fold(seq);
        count++;
    }
    
    return 0;
    
}

    
int
ViennaUnittest::run() {
    return 0;
}
    
int
ViennaUnittest::run_all() {
    String name = "ViennaUnittest";
    typedef int (ViennaUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_creation"        ] = &ViennaUnittest::test_creation;
    func_map["test_folding"         ] = &ViennaUnittest::test_folding;
    func_map["test_folding_no_reset"] = &ViennaUnittest::test_folding_no_reset;
    func_map["test_bp_probs"        ] = &ViennaUnittest::test_bp_probs;
    //func_map["test_memory_leak"       ] = &ViennaUnittest::test_memory_leak;
    
    
    int failed = 0;
    for(auto const & kv : func_map) {
        try {
            (this->*kv.second)();
        }
        catch(std::exception const & e) {
            std::cout << name << "::" << kv.first << " returned ERROR! : " << e.what() << std::endl;
            failed++;
        }
    }
    return failed;
}
    
}
}