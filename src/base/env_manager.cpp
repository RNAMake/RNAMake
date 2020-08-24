

#include<base/env_manager.h>
namespace base {

void
EnvManager::set_envs() {
    
    const auto cwd = std::__fs::filesystem::current_path().string();  
    
    
    for(const auto& env : env_vars_) {
        const auto current_var = getenv(env.c_str());

        if(current_var != nullptr) {
             

        } else {
            if(env.find("RNAMAKE") != std::string::npos 
                    || env.find("X3DNA") != std::string::npos) {

            } else {

            }

        }

    }

}

}

