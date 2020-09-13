//
// Created by Joseph Yesselman on 9/12/20.
//

#ifndef RNAMAKE_APPS_BUILD_MOTIF_GRAPH_BUILD_MOTIF_GRAPH_HXX
#define RNAMAKE_APPS_BUILD_MOTIF_GRAPH_BUILD_MOTIF_GRAPH_HXX

#include <CLI/CLI.hpp>

#include "base/application.hpp"
#include "base/log.h"
#include "resources/resource_manager.h"

class BuildMotifGraph : base::Application {
public:

  BuildMotifGraph();

public: // application functions

  void
  setup_options() override ;

  void
  parse_command_line(
    int,
    const char **) override ;

  void
  run() override ;

public: // getters

  base::LogLevel
  log_level() const {
    return base::log_level_from_str(_parameters.log_level);
  }

public:
  struct Parameters {
    String log_level = "info";
    String pdbs = "";
    String build_file = "";
  };

public:
  CLI::App _app;

private:
  // must be initialized a runtime
  resources::Manager & _rm;

  Parameters _parameters;

};

#endif //RNAMAKE_APPS_BUILD_MOTIF_GRAPH_BUILD_MOTIF_GRAPH_HXX
