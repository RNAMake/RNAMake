//
//  settings.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>

//RNAMake Headers
#include "base/types.h"
#include "base/settings.h"

String
base_dir() {

    char* base_path;
    base_path = std::getenv ("RNAMAKE");
    if (base_path==NULL) {
        std::cout << "cannot find environmental path RNAMAKE, please set it" << std::endl;
        exit(EXIT_FAILURE);
    }
    return String(base_path);
}

String
resources_path() {
    String base_path = base_dir();
    return base_path + "/rnamake/resources/";
}

String
lib_path() {
    String base_path = base_dir();
    return base_path + "/rnamake/lib/RNAMake/";
}

String
motif_dirs() {
    String base_path = base_dir();
    return base_path + "/rnamake/resources/motifs/";
}


String
x3dna_path() {
    String path = resources_path() + "x3dna/osx/";
    return path;
    
}

String
unittest_resource_dir() {
    return base_dir() + "/rnamake/lib/RNAMake/unittests_new/unittest_resources/";
}


