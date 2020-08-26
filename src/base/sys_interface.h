//
//  sys_interface.h
//  RNAMake
//
//  Created by Chris Jurich on Aug 11 2020
//  Copyright (c) 2020 Joseph Yesselman. All rights reserved.
//
#ifndef __SYS_INTERFACE_H__
#define __SYS_INTERFACE_H__

#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <cstdlib>


#include <nlohmann/json.hpp>
#include <base/types.h>

namespace base {

String execute_command( const char* );
String execute_command( String const &);

nlohmann::json
execute_command_json( const char* );

nlohmann::json
execute_command_json( String const &);

void
json_cleanup();

}
#endif // __SYS_INTERFACE_H__
