#include "../common.hpp"

#include <base/types.h>
#include <base/env_manager.h>

TEST_CASE( "Test Options for storing options for classes", "[Options]" ) {
    const auto vars = Strings{"RNAMAKE","X3DNA" };    
    auto EM = base::EnvManager{vars};
    
    EM.set_envs();
}
