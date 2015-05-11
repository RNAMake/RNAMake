//
//  x3dna_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "x3dna_unittest.h"
#include "util/file_io.h"


//get_ref_frame_path is private
/*int
X3dnaUnittest::test_get_ref_frame() {
    String m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6";
    
    String path = x_._get_ref_frame_path(m_path);
    std::cout << path << std::endl;
    return 1;
}*/

int
X3dnaUnittest::test_generate_ref_frame() {
    String m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6_2";
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
    String m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6_2";
    x_.generate_dssr_file(m_path);
    
    if(!(file_exists("p4p6_2_dssr.out"))) {
        throw "did not produce p4p6_2_dssr.out\n";
    }
    
    std::remove("p4p6_2_dssr.out");
    
    return 1;
}



int
X3dnaUnittest::run() {
    //if (test_get_ref_frame() == 0)    {  std::cout << "test_get_ref_frame failed" << std::endl; }
    if (test_generate_ref_frame() == 0)    {  std::cout << "test_generate_ref_frame failed" << std::endl; }
    if (test_generate_dssr_file() == 0)    {  std::cout << "test_generate_ref_frame failed" << std::endl; }

    return 1;
    
}