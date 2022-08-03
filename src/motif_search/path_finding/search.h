//
// Created by Joseph Yesselman on 3/17/19.
//

#ifndef RNAMAKE_NEW_PATH_FINDING_SEARCH_H
#define RNAMAKE_NEW_PATH_FINDING_SEARCH_H

#include <base/option.h>
#include <motif/motif_state_aligner.h>
#include <motif_search/path_finding/node.h>
#include <motif_search/path_finding/scorer.h>
#include <motif_search/path_finding/selector.h>
#include <motif_search/search.h>
#include <motif_search/solution_filter.h>

namespace motif_search {
namespace path_finding {

struct Parameters {
  bool sterics, helix_end;
  int max_node_level, min_size, max_size, max_solutions, min_node_level;
  float accept_score, min_ss_score, max_steps;
  bool return_best;
};

class Search : public motif_search::Search {
public:
  Search(ScorerOP scorer, SelectorOP selector, SolutionFilterOP solution_filter)
      : motif_search::Search("path_finding"), _scorer(scorer->clone()),
        _selector(selector->clone()),
        _solution_filter(solution_filter->clone()),
        _aligner(motif::MotifStateAligner()) {
    _motif_names = Strings();
    _motif_names.reserve(100);
    _parameters = Parameters();
    setup_options();
    update_var_options();
  }

  ~Search() {}

  motif_search::Search *clone() const { return new Search(*this); };

public:
  void setup(motif_search::ProblemOP);

  void start() {}

  bool finished() { return false; }

  motif_search::SolutionOP next();

protected:
  void update_var_options();

  void setup_options();

private:
  motif_data_structure::MotifStateGraphOP _graph_from_node(NodeOP);

  inline bool _accept_node(Node const &n) {
    if (n.ss_score() > _parameters.min_ss_score) {
      return false;
    }
    if (n.level() < _parameters.min_node_level) {
      return false;
    }
    // this is bad ... at motif_type to MotifState? -- JDY
    if (_parameters.helix_end && n.state()->name()[0] != 'H') {
      return false;
    }
    return true;
  }

  bool _steric_clash(motif::MotifState const &, Node const &);

  void _get_solution_motif_names(NodeOP, Strings &);

private:
  ScorerOP _scorer;
  SelectorOP _selector;
  SolutionFilterOP _solution_filter;
  Parameters _parameters;
  NodeQueue _queue;
  motif::MotifStateAligner _aligner;
  util::StericLookupNewOP _lookup;

  Strings _motif_names;
  bool _using_lookup, _enumerating;
};

} // namespace path_finding
} // namespace motif_search

#endif // RNAMAKE_NEW_PATH_FINDING_SEARCH_H
