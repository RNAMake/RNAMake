

#include<base/env_manager.h>

namespace base {

void
EnvManager::set_envs() {
    char buffer[500];
    if(getcwd(buffer,500) == nullptr) {
        std::cout<<"Error in getting current working directory\nExiting\n";
        exit(1);
    }
    const auto cwd = String{buffer};  
    
    
    for(const auto& env : env_vars_) {
        const auto current_var = getenv(env.c_str());
        if(current_var != nullptr) { continue; }
        
        if(env.find("RNAMAKE") != std::string::npos ) {
                const auto rnamake_it = cwd.find("RNAMake") ;
                std::cout<<"\"RNAMAKE\" environment variable not set...\n";
            if(rnamake_it != std::string::npos) {
                const auto rnamake_var = cwd.substr(0,rnamake_it + 7);
                std::cout<<"\tattempting to set equal to \""<<rnamake_var<<"\"...\n";
                const auto error_code = setenv("RNAMAKE",rnamake_var.c_str(),1);     
                
                if(error_code != 0) {
                    std::cout<<"\tAttempt to update \"RNAMAKE\" failed.\n\tPlease set variable and re-run application. Exiting.\n";
                    exit(1);
                } else {
                    std::cout<<"\tUpdate successful! Continuing...\n";
                }
                
            } else {
                std::cout<<"\tCannot determine \"RNAMAKE\" variable automatically.\n\tPlease set this variable and re-run the application. Exiting.\n";
                exit(1);
            }

        } else if ( env.find("X3DNA") != std::string::npos) {
                const auto rnamake_it = cwd.find("RNAMake") ;
                std::cout<<"\"X3DNA\" environment variable not set...\n";
            if(rnamake_it != std::string::npos) {
                //const auto x3dna_var = cwd.substr(0,rnamake_it + 7) + "/resources/" + base::get_os ;
                const auto x3dna_var = base::x3dna_path(); 
                std::cout<<"\tattempting to set equal to \""<<x3dna_var<<"\"...\n";
                const auto error_code = setenv("X3DNA",x3dna_var.c_str(),1);     
                
                if(error_code != 0) {
                    std::cout<<"\tAttempt to update \"X3DNA\" failed.\n\tPlease set variable and re-run application. Exiting.\n";
                    exit(1);
                } else {
                    std::cout<<"\tUpdate successful! Continuing...\n";
                }
                
            } else {
                std::cout<<"\tCannot determine \"X3DNA\" variable automatically.\n\tPlease set this variable and re-run the application. Exiting.\n";
                exit(1);
            }


        } else {

        }

        

    }

}

}

