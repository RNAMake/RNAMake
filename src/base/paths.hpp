//
//  settings.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__settings__
#define __RNAMake__settings__

#include <cstdio>

// RNAMake Headers
#include <base/types.hpp>

namespace base::path {

// This all needs to be replaced by
// https://en.cppreference.com/w/cpp/filesystem/path

String filename(String const &);

String parent_dir(String const &);

String rnamake_path();

String resources_path();

String unittest_resource_path();

void get_lines_from_file(String const &,
                         Strings & /* return */) noexcept(false);

} // namespace base::path

#endif /* defined(__RNAMake__settings__) */
