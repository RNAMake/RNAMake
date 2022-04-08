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

/* TODO implement these in paths with new
String filename(String const &);

String base_dir(String const &);
*/

String rnamake_path();

String resources_path();

String unittest_resource_path();

} 

#endif /* defined(__RNAMake__settings__) */
