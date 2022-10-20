//
// Created by Joseph Yesselman on 11/4/17.
//

#ifndef TEST_MOTIF_STATE_MONTE_CARLO_H
#define TEST_MOTIF_STATE_MONTE_CARLO_H

#include <base/option.h>
#include <motif_data_structure/motif_state_graph.hpp>
#include <motif_data_structure/motif_state_tree.h>
#include <util/monte_carlo.h>
#include <util/steric_lookup.hpp>

namespace motif_search {

struct MotifStateMonteCarloSolution {
  inline MotifStateMonteCarloSolution(motif_data_structure::MotifGraphOP n_mg,
                                      float n_score)
      : mg(n_mg), score(n_score) {}

  motif_data_structure::MotifGraphOP mg;
  float score;
};

struct MotifStateMonteCarloSolutionNew {
  inline MotifStateMonteCarloSolutionNew(
      motif_data_structure::MotifStateGraphOP n_msg, float n_score)
      : msg(n_msg), score(n_score) {}

  motif_data_structure::MotifStateGraphOP msg;
  float score;
};

typedef std::shared_ptr<MotifStateMonteCarloSolution>
    MotifStateMonteCarloSolutionOP;
typedef std::shared_ptr<MotifStateMonteCarloSolutionNew>
    MotifStateMonteCarloSolutionNewOP;

class MotifStateMonteCarlo {
public:
  MotifStateMonteCarlo(std::vector<motif::MotifStateOPs> const &mses)
      : mses_(mses), mc_(util::MonteCarlo(0.5f)),
        rng_(util::RandomNumberGenerator()), options_(base::Options()) {
    using_lookup_ = 0;
    setup_options();
  }

public:
  void setup(motif_data_structure::MotifStateGraphOP msg, int, int, int, int,
             bool);

  void run();

  void start();

  MotifStateMonteCarloSolutionOP next();

  MotifStateMonteCarloSolutionNewOP next_state();

  bool finished();

public:
  inline void lookup(util::StericLookupNew const &sl) {
    using_lookup_ = 1;
    lookup_ = sl;
  }

private:
  double get_score(structure::BasepairStateOP);

  float perform_motif_swap(float);

  void perform_motif_swap_no_clash();

  bool _steric_clash(motif_data_structure::MotifStateGraphOP);

  bool _seen_solution(motif_data_structure::MotifStateGraphOP);

protected:
  void setup_options();

  void update_var_options();

public: // option wrappers
  inline base::Options &options() { return _options; }

  inline float get_int_option(String const &name) {
    return _options.get_int(name);
  }

  inline float get_float_option(String const &name) {
    return _options.get_float(name);
  }

  inline String get_string_option(String const &name) {
    return _options.get_string(name);
  }

  inline bool get_bool_option(String const &name) {
    return _options.get_bool(name);
  }

  inline bool has_option(String const &name) {
    return _options.has_option(name);
  }

  template <typename T>
  void set_option_value(String const &name, T const &val) {
    _options.set_value(name, val);
    update_var_options();
  }

private:
  util::MonteCarlo _mc;
  util::RandomNumberGenerator _rng;
  util::StericLookup _lookup;
  math::Vector3s _beads;
  base::Options _options;
  std::vector<motif::MotifStateOPs> _mses;
  structure::BasepairStateOP _end, _end_flip;
  motif::MotifStateOP _start_m;
  motif_data_structure::MotifStateGraphOP _msg;
  double _score, _r_diff, _r_diff_flip;
  int _ni, _nj, _ei, _ej;
  int _org_num;
  int _stage;
  int _step;
  bool _target_an_aligned_end;
  StringIntMap _seen;
  bool _enumerating;
  int _using_lookup;
  // options
  bool _sterics;
  int _max_solutions;
  int _stages;
  int _steps;
  float _accept_score;
};

typedef std::shared_ptr<MotifStateMonteCarlo> MotifStateMonteCarloOP;

} // namespace motif_search

#endif // TEST_MOTIF_STATE_MONTE_CARLO_H
