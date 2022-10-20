//
// Created by Joseph Yesselman on 2019-03-30.
//

#ifndef RNAMAKE_NEW_EXHAUSTIVE_SEARCH_H
#define RNAMAKE_NEW_EXHAUSTIVE_SEARCH_H

#include <motif_search/exhaustive/motif_state_enumerator.h>
#include <motif_search/exhaustive/scorer.h>
#include <motif_search/search.h>
#include <motif_search/solution_filter.h>
#include <motif_search/solution_topology.h>

namespace motif_search {
namespace exhaustive {

class Search : public motif_search::Search {
public:
  struct Parameters {
    bool sterics, helix_end;
    int max_node_level, min_size, max_size, max_solutions, min_node_level;
    float accept_score, min_ss_score, max_steps;
    bool return_best;
  };

public:
  Search(ScorerOP scorer, SolutionToplogy const &sol_top,
         SolutionFilterOP filter)
      : motif_search::Search("exhaustive"), _scorer(scorer->clone()),
        _enumerator(MotifStateEnumerator(sol_top)), _filter(filter->clone()),
        _parameters(Parameters()) {
    _finished = false;
    setup_options();
    update_var_options();
  }

  ~Search() {}

  motif_search::Search *clone() const { return new Search(*this); };

public:
  void setup(ProblemOP p) {
    _enumerator.set_size_limit(get_int_option("max_size"));
    _enumerator.start(p->start);
    _scorer->set_target(p->end, p->target_an_aligned_end);
  }

  void start() {}

  bool finished() { return _finished; }

  SolutionOP next() {
    auto best = 1000.0f;
    while (!_enumerator.finished()) {
      _current = _enumerator.top_state();
      _score = _scorer->score(*_scorer->end_states()[1]);

      if (_score < _parameters.accept_score) {
        auto msg = _get_graph_from_solution();
        _enumerator.next();
        _get_solution_motif_names(msg);
        if (!_filter->accept(_motif_names)) {
          continue;
        }

        return std::make_shared<motif_search::Solution>(msg, _score);
      }
      _enumerator.next();
    }

    _finished = true;
    return SolutionOP(nullptr);
  }

private:
  void _get_solution_motif_names(motif_data_structure::MotifStateGraphOP msg) {
    _motif_names.resize(0);
    for (auto const &n : *msg) {
      _motif_names.push_back(n->data()->name());
    }
  }

  motif_data_structure::MotifStateGraphOP _get_graph_from_solution();

protected:
  void setup_options();

  void update_var_options();

private:
  Parameters _parameters;
  ScorerOP _scorer;
  SolutionFilterOP _filter;
  MotifStateEnumerator _enumerator;
  motif::MotifStateOP _current;
  float _score;
  bool _finished;
  Strings _motif_names;
};

} // namespace exhaustive
} // namespace motif_search

#endif // RNAMAKE_NEW_EXHAUSTIVE_SEARCH_H
