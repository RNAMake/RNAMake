//
//  x3dna_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "x3dna_unittest.h"
#include "util/file_io.h"


int
X3dnaUnittest::test_generate_ref_frame() {
    String m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6_2/p4p6_2.pdb";
    x_.generate_ref_frame(m_path);
    
    if(!(file_exists("ref_frames.dat"))) {
        throw "did not produce ref_frames.dat\n";
    }
       
    if(file_exists("basepairs.pdb")) {
        throw "produced basepairs.pdb which is not needed\n";
    }
    
    //cleanup
    std::remove("ref_frames.dat");
    
    return 1;
}

int
X3dnaUnittest::test_generate_dssr_file() {
    String m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6_2/p4p6_2.pdb";
    x_.generate_dssr_file(m_path);
    
    if(!(file_exists("p4p6_2_dssr.out"))) {
        throw "did not produce p4p6_2_dssr.out\n";
    }
    
    std::remove("p4p6_2_dssr.out");
    
    return 1;
}

int
X3dnaUnittest::test_res_compare() {
    X3Residue res1(1, "A", "");
    X3Residue res2(1, "A", "");
    X3Residue res3(2, "A", "");
    X3Residue res4(1, "B", "");

    if(!(res1 == res2)) { return 0; }
    if(res1 == res3)    { return 0; }
    if(res1 == res4)    { return 0; }

    return 1;
}

int
X3dnaUnittest::test_get_basepairs() {
    String m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6/p4p6.pdb";
    auto basepairs = x_.get_basepairs(m_path);
    
    if(file_exists("p4p6_dssr.out") || file_exists("ref_frames.dat")) {
        throw UnittestException("produced x3dna files when they already exists");
    }
    
    m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6_2/p4p6_2.pdb";
    basepairs = x_.get_basepairs(m_path);

    if(! file_exists("p4p6_2_dssr.out") || ! file_exists("ref_frames.dat")) {
        throw UnittestException("did not produce x3dna files when required");
    }
    
    //cleanup
    std::remove("ref_frames.dat");
    std::remove("p4p6_2_dssr.out");

    return 1;
}

int
X3dnaUnittest::test_get_motifs() {
    String m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6/p4p6.pdb";
    auto motifs = x_.get_motifs(m_path);
    
    
    return 1;
}



int
X3dnaUnittest::run() {
    //if (test_generate_ref_frame() == 0)    {  std::cout << "test_generate_ref_frame failed" << std::endl; }
    //if (test_generate_dssr_file() == 0)    {  std::cout << "test_generate_ref_frame failed" << std::endl; }
    //if (test_res_compare() == 0)           {  std::cout << "test_res_compare failed" << std::endl;}
    //if (test_get_basepairs() == 0)         {  std::cout << "test_get_basepairs failed" << std::endl;}
    if (test_get_motifs() == 0)            {  std::cout << "test_get_motifs failed" << std::endl;}

    return 1;
    
}


int
X3dnaUnittest::run_all() {
    String name = "X3dnaUnittest";
    typedef int (X3dnaUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_generate_ref_frame" ] = &X3dnaUnittest::test_generate_ref_frame;
    func_map["test_generate_dssr_file" ] = &X3dnaUnittest::test_generate_dssr_file;
    func_map["test_res_compare"        ] = &X3dnaUnittest::test_res_compare;
    func_map["test_get_basepairs"      ] = &X3dnaUnittest::test_get_basepairs;
    
    for(auto const & kv : func_map) {
        try {
            int result = (this->*kv.second)();
            if(result == 0) {
                std::cout << name << "::" << kv.first << " FAILED!" << std::endl;
            }
        }
        catch(...) {
            std::cout << name << "::" << kv.first << " returned ERROR!" << std::endl;
        }
        
    }
    
    return 0;
}





