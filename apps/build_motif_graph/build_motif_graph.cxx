//
// Created by Joseph Yesselman on 9/12/20.
//

#include <cstdlib>

// app header
#include <thermo_fluctuation/graph/simulation.h>

#include "build_motif_graph/build_motif_graph.hxx"

// rnamake headers
#include "base/backtrace.h"
#include "util/csv.h"

// init
// ///////////////////////////////////////////////////////////////////////////////////////////

BuildMotifGraph::BuildMotifGraph()
    : base::Application(), _rm(resources::Manager::instance()) {}

// app functions
// //////////////////////////////////////////////////////////////////////////////////

void BuildMotifGraph::setup_options() {
  _app.add_option_group("Core Inputs");
  _app.add_option("--build", _parameters.build_file,
                  "file with instructions for how to build graph")
      ->required()
      ->group("Core Inputs");

  _app.add_option("--pdbs", _parameters.pdbs,
                  "path to a PDB file with input RNA structure")
      ->default_val("")
      ->group("Core Inputs");

  _app.add_option("--ensembles", _parameters.ensemble_dirs, "path to ensemble ")
      ->default_val("")
      ->group("Core Inputs");

  _app.add_option("--connect", _parameters.connect,
                  "file with instructions for how to build graph")
      ->group("Core Inputs");

  _app.add_option("--seq", _parameters.sequence, "")->group("Core Inputs");
  _app.add_option("--scorer", _parameters.scorer, "");
  _app.add_option("--steps", _parameters.steps, "");
  _app.add_option("--runs", _parameters.runs, "");
  _app.add_option("--resample-sequence", _parameters.resample_sequence);

  _app.add_flag("--output_pdb", _parameters.output_pdb);
  _app.add_flag("--sterics", _parameters.sterics);
}

void BuildMotifGraph::parse_command_line(int argc, const char** argv) {
  base::init_logging();
  // handle additional pdb files supplied
  _parse_pdb_files();

  // handle additional ensembles supplied
  _parse_ensemble_files();
}

void BuildMotifGraph::run() {
  _mg = motif_data_structure::MotifGraph();
  _build_motif_graph_from_csv(_mg);

  if (_parameters.connect.empty()) {
    LOG_INFO << "completed, writing pdb to \"test.pdb\"";
    std::cout << _mg.sequence() << std::endl;
    std::cout << _mg.dot_bracket() << std::endl;
    _mg.write_pdbs();
    _mg.to_pdb("test.pdb", 1, 1);
    exit(0);
  }
  auto spl = base::split_str_by_delimiter(_parameters.connect, ",");
  auto ni = std::stoi(spl[0]);
  auto end_name = spl[1];
  auto nie = _mg.get_node(ni)->data()->get_end_index(end_name);
  _mg.add_connection(_mg.last_node()->index(), std::stoi(spl[0]),
                    _mg.last_node()->data()->end_name(1), spl[1]);

  if (_parameters.sequence.empty() && _parameters.resample_sequence == -1) {
    LOG_INFO << "completed, writing pdb to \"test.pdb\"";
    std::cout << _mg.sequence() << std::endl;
    std::cout << _mg.dot_bracket() << std::endl;
    _mg.to_pdb("test.pdb", 1, 1);
    exit(0);
  }

  _mg.replace_ideal_helices();
  _run_thermo_fluc(_mg);


  // best_mg->write_pdbs();
  // mg.to_pdb("test.pdb", 1, 1);
}

// app private functions
// ///////////////////////////////////////////////////////////////////////////

void BuildMotifGraph::_parse_pdb_files() {
  auto pdb_files = base::split_str_by_delimiter(_parameters.pdbs, ",");
  if (!pdb_files.empty()) {
    LOG_INFO << "loading pdbs from --extra_pdb flag: " << _parameters.pdbs;
  }
  for (auto const& pdb_file : pdb_files) {
    if (!base::file_exists(pdb_file)) {
      LOG_ERROR << "invalid pdb path: " << pdb_file;
      exit(1);
    }
    LOG_INFO << "loading " << pdb_file;
    _rm.add_motif(pdb_file, "", util::MotifType::TCONTACT);
  }
}

void BuildMotifGraph::_parse_ensemble_files() {
  auto ensemble_files =
      base::split_str_by_delimiter(_parameters.ensemble_dirs, ",");
  if (!ensemble_files.empty()) {
    LOG_INFO << "loading ensemble files from --ensembles "
             << _parameters.ensemble_dirs;
  }
  for (auto const& ensemble_file : ensemble_files) {
    if (!base::file_exists(ensemble_file)) {
      LOG_ERROR << "invalid ensemble path: " << ensemble_file;
      exit(1);
    }
    auto mes = std::vector<motif::MotifEnsembleOP>();
    motif::motif_ensemble_from_csv_file(ensemble_file, mes);
    for (auto const& me : mes) {
      auto name = me->members()[0]->motif->name();
      auto end_name = me->members()[0]->motif->end_name(0);
      _rm.register_motif_ensemble(name, end_name, me);
    }
  }
}

void BuildMotifGraph::_build_motif_graph_from_csv(
    motif_data_structure::MotifGraph& mg) {
  auto in = io::CSVReader<5>(_parameters.build_file);
  mg.set_option_value("sterics", false);
  String motif, align_end, parent_end_name;
  int parent, parent_end_index;
  in.read_header(io::ignore_missing_column, "motif", "align_end", "parent",
                 "parent_end_index", "parent_end_name");
  while (in.read_row(motif, align_end, parent, parent_end_index,
                     parent_end_name)) {
    auto m = motif::MotifOP(nullptr);
    if (align_end.empty()) {
      m = _rm.motif(motif);
    } else {
      m = _rm.motif(motif, "", align_end);
    }
    auto spl = base::split_str_by_delimiter(m->name(), ".");
    if (spl[0] == "HELIX") {
      m->mtype(util::MotifType::HELIX);
    }

    if (!parent_end_name.empty()) {
      mg.add_motif(m, parent, parent_end_name);
    } else {
      mg.add_motif(m);
    }
  }
}

void BuildMotifGraph::_run_thermo_fluc(
    motif_data_structure::MotifGraph & mg) {
  // save some characters with namespace
  namespace tfg = thermo_fluctuation::graph;
  auto spl = base::split_str_by_delimiter(_parameters.connect, ",");
  auto ni = std::stoi(spl[0]);
  auto end_name = spl[1];
  auto nie = _mg.get_node(ni)->data()->get_end_index(end_name);
  auto last_m = _mg.get_node(ni)->connections()[nie]->partner(ni);
  _mg.replace_helical_sequence(_parameters.sequence);
  auto index_hash = std::map<int, int>();
  auto mseg = motif_data_structure::MotifStateEnsembleGraph();
  motif_data_structure::motif_state_ensemble_graph_from_motif_graph(
      mg, _rm, mseg, index_hash);
  auto start = data_structure::NodeIndexandEdge{index_hash[ni], nie};
  auto end = data_structure::NodeIndexandEdge{index_hash[last_m->index()], 1};

  auto thermo_scorer = tfg::ScorerOP();
  auto sterics = tfg::sterics::StericsOP();
  if (_parameters.scorer == "OldFrameScorer") {
    thermo_scorer = tfg::ScorerOP(std::make_shared<tfg::OldFrameScorer>());
  } else if (_parameters.scorer == "FrameScorer") {
    thermo_scorer = tfg::ScorerOP(std::make_shared<tfg::FrameScorer>());
  } else {
    LOG_ERROR << "invalid scorer: " << _parameters.scorer;
    exit(1);
  }

  if (_parameters.sterics) {
    sterics = tfg::sterics::StericsOP(
        std::make_shared<tfg::sterics::SelectiveSterics>(
            Ints{0}, Ints{index_hash[last_m->index()]}, 2.2f));
  } else {
    sterics =
        tfg::sterics::StericsOP(std::make_shared<tfg::sterics::NoSterics>());
  }
  auto avg = 0;
  auto best_mg = motif_data_structure::MotifGraphOP();
  for (int i = 0; i < _parameters.runs; i++) {
    auto thermo_sim_ =
        std::make_shared<tfg::Simulation>(thermo_scorer, sterics);
    thermo_sim_->set_option_value("cutoff", _parameters.cutoff);
    thermo_sim_->setup(mseg, start, end);
    thermo_sim_->next();
    if (i == 0) {
      best_mg = thermo_sim_->get_motif_graph();
    }
    auto count = 0;
    auto best = thermo_sim_->get_score();
    std::cout << best << std::endl;
    auto under_cutoff = 0;

    for (int s = 0; s < _parameters.steps; s++) {
      under_cutoff = thermo_sim_->next();
      if (under_cutoff) {
        count += 1;
      }
      if (thermo_sim_->get_score() < best && _parameters.output_pdb) {
        // LOG_INFO << "best dist score: " << best;
        best = thermo_sim_->get_score();
        if (_parameters.output_pdb) {
          best_mg = thermo_sim_->get_motif_graph();
        }
      }
    }
    avg += count;
  }
  std::cout << avg / _parameters.runs << std::endl;
  // connection is not perserved through simulation ... bring it back
  if (_parameters.output_pdb) {
    best_mg->add_connection(start.node_index, end.node_index,
                            best_mg->get_node(start.node_index)
                                ->data()
                                ->ends()[start.edge_index]
                                ->name(),
                            best_mg->get_node(end.node_index)
                                ->data()
                                ->ends()[end.edge_index]
                                ->name());
    best_mg->to_pdb("test.pdb", 1, 1);
  }
}

int main(int argc, const char** argv) {
  // must add this for all apps!
  std::set_terminate(base::print_backtrace);

  auto app = BuildMotifGraph();
  app.setup_options();
  CLI11_PARSE(app._app, argc, argv);
  app.parse_command_line(argc, argv);
  app.run();
  return 0;
}
