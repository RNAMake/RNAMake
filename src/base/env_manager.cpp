

#include <base/env_manager.h>

#include "base/settings.h"

namespace base {

void EnvManager::set_envs() {
  char buffer[500];
  if (getcwd(buffer, 500) == nullptr) {
    LOG_ERROR << "Error in getting current working directory\nExiting";
    exit(1);
  }
  const auto cwd = String{buffer};

  for (const auto& env : env_vars_) {
    const auto current_var = getenv(env.c_str());
    if (current_var != nullptr) {
      continue;
    }

    if (env.find("RNAMAKE") != std::string::npos) {
      const auto rnamake_it = cwd.find("RNAMake");
      LOG_WARNING << "\"RNAMAKE\" environment variable not set...exiting";
      exit(1);
      /*if(rnamake_it != std::string::npos) {
          const auto rnamake_var = cwd.substr(0,rnamake_it + 7);
          LOG_WARNING<<"\tattempting to set equal to \""<<rnamake_var<<"\"...";
          const auto error_code = setenv("RNAMAKE",rnamake_var.c_str(),1);

          if(error_code != 0) {
              LOG_ERROR<<"\tAttempt to update \"RNAMAKE\" failed.\n\tPlease set
      variable and re-run application. Exiting."; exit(1); } else {
              LOG_WARNING<<"\tUpdate successful! Continuing...";
          }

      } else {
          std::cout<<"\tCannot determine \"RNAMAKE\" variable
      automatically.\n\tPlease set this variable and re-run the application.
      Exiting."; exit(1);
      } */

    } else if (env.find("X3DNA") != std::string::npos) {
      const auto rnamake_it = cwd.find("RNAMake");
      LOG_WARNING << "\"X3DNA\" environment variable not set... exiting";
      exit(0);
      /*if(rnamake_it != std::string::npos) {
          const auto x3dna_var = base::x3dna_path();
          LOG_WARNING<<"\tattempting to set equal to \""<<x3dna_var<<"\"...";
          const auto error_code = setenv("X3DNA",x3dna_var.c_str(),1);

          if(error_code != 0) {
              LOG_ERROR<<"\tAttempt to update \"X3DNA\" failed.\n\tPlease set
      variable and re-run application. Exiting."; exit(1); } else {
              LOG_WARNING<<"\tUpdate successful! Continuing...";
          }

      } else {
          LOG_ERROR<<"\tCannot determine \"X3DNA\" variable
      automatically.\n\tPlease set this variable and re-run the application.
      Exiting."; exit(1);
      }
      */

    } else {
    }
  }
}

}  // namespace base
