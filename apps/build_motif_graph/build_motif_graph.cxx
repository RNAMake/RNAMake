//
// Created by Joseph Yesselman on 9/12/20.
//

// app header
#include "build_motif_graph/build_motif_graph.hxx"

// rnamake headers
#include "base/backtrace.h"
#include "util/csv.h"
#include "motif_data_structure/motif_graph.h"

// init  ///////////////////////////////////////////////////////////////////////////////////////////

BuildMotifGraph::BuildMotifGraph():
  base::Application(),
  _rm(resources::Manager::instance()){
}

// app functions  //////////////////////////////////////////////////////////////////////////////////

void
BuildMotifGraph::setup_options() {
  _app.add_option_group("Core Inputs");
  _app.add_option("--pdbs", _parameters.pdbs, "path to a PDB file with input RNA structure")
    ->default_val("")
    ->group("Core Inputs");

  _app.add_option(
      "--build", _parameters.build_file, "file with instructions for how to build graph")
    ->required()
    ->group("Core Inputs");

}

void
BuildMotifGraph::parse_command_line(
  int argc,
  const char**argv) {
}

void
BuildMotifGraph::run() {
  auto pdb_files = base::split_str_by_delimiter(_parameters.pdbs, ",");
  if(!pdb_files.empty()) {
    LOG_INFO << " loading pdbs from --extra_pdb flag: " << _parameters.pdbs;
  }
  for(auto const & pdb_file : pdb_files) {
    if(!base::file_exists(pdb_file)) {
      LOG_ERROR << "invalid pdb path: " << pdb_file;
      exit(0);
    }
    LOG_INFO << "loading " << pdb_file;
    _rm.add_motif(pdb_file);
  }

  auto csv = util::csv::Reader();
  auto csv_table = csv.read_csv(_parameters.build_file);
  auto mg = motif_data_structure::MotifGraph();

  for(auto const & row : *csv_table) {
    std::cout << row.get_int_val("motif_name") << std::endl;
  }

}

int
main(
  int argc,
  const char**argv) {
  //must add this for all apps!
  std::set_terminate(base::save_backtrace);

  auto app = BuildMotifGraph();
  app.setup_options();
  CLI11_PARSE(app._app, argc, argv);
  //start logging
  base::init_logging(app.log_level());
  app.run();
  return 0;
}
