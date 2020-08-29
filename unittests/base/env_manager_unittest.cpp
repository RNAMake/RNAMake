#include "../common.hpp"

#include <cstdlib>

#include <base/types.h>
#include <base/log.h>
#include <base/env_manager.h>

TEST_CASE( "Test environment variable setting", "[EnvManager]" ) {
    
    base::init_logging(base::LogLevel::FATAL);

    auto unset_vars(0);
    const auto rnamake_orig = std::getenv("RNAMAKE");
    
    if(rnamake_orig != nullptr) {
        LOGW<<"NOTE: Unsetting environment variable \"RNAMAKE\" for unittesting ...";
        unsetenv("RNAMAKE");
        auto rnamake_setter = base::EnvManager{Strings{"RNAMAKE"}};
        rnamake_setter.set_envs();
        const auto rnamake_updated = std::getenv("RNAMAKE");

        LOGW<<"NOTE: Environment variable \"RNAMAKE\" has been restored";
        
        REQUIRE(String{rnamake_orig} == String{rnamake_updated});

    } else {
        ++unset_vars;
    }

    const auto x3dna_orig = std::getenv("X3DNA");

    if(x3dna_orig != nullptr) {
        
        LOGW<<"NOTE: Unsetting environment variable \"X3DNA\" for unittesting ...";
        unsetenv("X3DNA");
        auto x3dna_setter = base::EnvManager{Strings{"X3DNA"}};
        x3dna_setter.set_envs();
        const auto x3dna_updated = std::getenv("X3DNA");
        
        LOGW<<"NOTE: Environment variable \"X3DNA\" has been restored";
        
        REQUIRE(String{x3dna_orig} == String{x3dna_updated});
    } else {
        ++unset_vars;
    }
    if(unset_vars) {
        LOGW<<"========================================================================"; 
        LOGW<<"WARNING: You must set the envionrment variables \"X3DNA\" and \"RNAMAKE\"";
        LOGW<<"To check if they are set, run:";
        LOGW<<"\t>echo $RNAMAKE";
        LOGW<<"\t\tor";
        LOGW<<"\t>echo $X3DNA";
        LOGW<<"in the terminal to see their current value(s).";
        LOGW<<"======================================================================="; 
    }
    
    REQUIRE(unset_vars == 0);    
}
