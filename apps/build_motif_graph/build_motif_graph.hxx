//
// Created by Joseph Yesselman on 9/12/20.
//

#ifndef RNAMAKE_APPS_BUILD_MOTIF_GRAPH_BUILD_MOTIF_GRAPH_HXX
#define RNAMAKE_APPS_BUILD_MOTIF_GRAPH_BUILD_MOTIF_GRAPH_HXX

#include <CLI/CLI.hpp>

#include "base/application.hpp"
#include "base/log.h"
#include "resources/resource_manager.h"
#include "motif_data_structure/motif_graph.h"
#include <thermo_fluctuation/graph/simulation.h>


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

private:
  void
  _parse_pdb_files();

  void
  _parse_ensemble_files();

  void
  build_motif_graph_from_csv(
    motif_data_structure::MotifGraph & /* return */);


public:
  struct Parameters {
    String log_level = "info";
    String pdbs = "";
    String build_file = "";
    String connect = "";
    String sequence = "";
    String ensemble_dirs = "";
    bool output_pdb = false;
  };

public:
  CLI::App _app;

private:
  // must be initialized a runtime
  resources::Manager & _rm;
  // default initiation
  Parameters _parameters = Parameters();

};

#endif //RNAMAKE_APPS_BUILD_MOTIF_GRAPH_BUILD_MOTIF_GRAPH_HXX
