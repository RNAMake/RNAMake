//
//  segmenter_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "segmenter_unittest.h"
#include "motif/pose_factory.h"
#include "motif_tools/segmenter.h"

int
SegmenterUnittest::test_creation() {
    /*Segmenter s;
    PoseFactory pf;
    String path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6/p4p6.pdb";
    auto p = pf.pose_from_file(path);
    auto bp1 = p->get_basepair("A112-A208")[0];
    auto bp2 = p->get_basepair("A132-A191")[0];
    s.apply(p, BasepairOPs { bp1, bp2} );*/
    return 1;
}

int
SegmenterUnittest::run() {
    if (test_creation() == 0)        { std::cout << "test_creation failed" << std::endl;  }
    return 1;
}