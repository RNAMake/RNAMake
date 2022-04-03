//
//  settings.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__settings__
#define __RNAMake__settings__

#include <stdio.h>

// RNAMake Headers
#include <base/types.hpp>

namespace base {

String get_os_name();

String base_dir();

String resources_path();

String lib_path();

String motif_dirs();

String unittest_resource_dir();

}  // namespace base

#endif /* defined(__RNAMake__settings__) */
