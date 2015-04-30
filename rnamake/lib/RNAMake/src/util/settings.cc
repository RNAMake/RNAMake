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
#include "util/settings.h"

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