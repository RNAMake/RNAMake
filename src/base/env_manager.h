//
//  sys_interface.h
//  RNAMake
//
//  Created by Chris Jurich on Aug 23 2020
//  Copyright (c) 2020 Joseph Yesselman. All rights reserved.
//

#ifndef __ENV_MANAGER_H__
#define __ENV_MANAGER_H__

#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <filesystem>

#include <base/types.h>

namespace base {

class EnvManager{
    private:
        Strings env_vars_;
    public:
        EnvManager(Strings const& env_vars ) : 
            env_vars_(env_vars) {

            }
    public: 
        void
        add_env(String const& env) {
            env_vars_.push_back(std::move(env));
        }

    public:
        void
        set_envs();
};

}



#endif // __ENV_MANAGER_H__
