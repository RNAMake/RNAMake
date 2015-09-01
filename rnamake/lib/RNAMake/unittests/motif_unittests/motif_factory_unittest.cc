//
//  motif_factory_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/9/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_factory_unittest.h"
#include "util/settings.h"
#include "motif/motif_factory.h"

int
MotifFactoryUnittest::test_creation() {
    MotifFactory mf;
    return 1;
}

int
MotifFactoryUnittest::test_load() {
    MotifFactory mf;
    String path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6/p4p6.pdb";
    mf.motif_from_file(path);
    return 1;
}

int
MotifFactoryUnittest::test_standardize_motif() {
    MotifFactory mf;
    auto path = motif_dirs() + "helices/HELIX.IDEAL";
    auto m = mf.motif_from_file(path);
    mf.standardize_motif(m);
    
    return 1;
}

int
MotifFactoryUnittest::test_can_align_motif_to_end() {
    MotifFactory mf;
    auto path = motif_dirs() + "helices/HELIX.IDEAL";
    auto m = mf.motif_from_file(path);
    auto m_added = mf.can_align_motif_to_end(m, 0);
    return 1;
}


int
MotifFactoryUnittest::run() {
    if (test_creation() == 0)            { std::cout << "test_creation failed" << std::endl;  }
    if (test_load() == 0)                { std::cout << "test_load failed" << std::endl;  }
    if (test_standardize_motif() == 0)   { std::cout << "test_standardize_motif failed" << std::endl;  }
    if (test_can_align_motif_to_end() == 0)   { std::cout << "test_can_align_motif_to_end failed" << std::endl;  }

    return 0;
}
