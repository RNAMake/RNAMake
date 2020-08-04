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
#include "base/log.h"

namespace base {

String
get_os_name() {
#ifdef _WIN32 || _WIN64
    return  String("Windows");
#elif __unix || __unix__
    return  String("unix");
#elif __APPLE__ || __MACH__
    return String("OSX");
#elif __linux__
    return String("Linux");
#endif

    throw std::runtime_error("cannot determine operating system");
}

String
base_dir() {

    char *base_path;
    base_path = std::getenv("RNAMAKE");
    if (base_path == NULL) {
        LOG_ERROR << "cannot find environmental path RNAMAKE, please set it" << std::endl;
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
    auto os_name = get_os_name();
    if (os_name == "OSX") { return resources_path() + "x3dna/osx/"; }
    if (os_name == "Linux" || os_name == "unix") { return resources_path() + "x3dna/linux/"; }
    throw std::runtime_error("unsupported operating system!");

}

String
unittest_resource_dir() {
    return base_dir() + "/rnamake/lib/RNAMake/unittests/unittest_resources/";
}

}


