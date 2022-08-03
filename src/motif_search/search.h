//
// Created by Joseph Yesselman on 3/17/19.
//

#ifndef RNAMAKE_NEW_SEARCH_H
#define RNAMAKE_NEW_SEARCH_H

#include <motif_data_structure/motif_state_graph.hpp>
#include <motif_search/problem.h>

namespace motif_search {

struct Solution {
  inline Solution(motif_data_structure::MotifStateGraphOP n_graph,
                  float n_score)
      : graph(n_graph), score(n_score) {}

  motif_data_structure::MotifStateGraphOP graph;
  float score;
};

typedef std::shared_ptr<Solution> SolutionOP;

class Search {
public:
public:
  Search(String const &name) : _name(name) {}

  virtual ~Search() = default;

  virtual Search *clone() const = 0;

public:
  virtual void setup(ProblemOP) = 0;

  virtual void start() = 0;

  virtual bool finished() = 0;

  virtual SolutionOP next() = 0;

public:
  String const &name() { return _name; }

public: // option wrappers
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

  template <typename T>
  void set_option_value(String const &name, T const &val) {
    _options.set_value(name, val);
    update_var_options();
  }

protected:
  virtual void setup_options() = 0;

  virtual void update_var_options() = 0;

protected:
  base::Options _options;
  String _name;
};

typedef std::shared_ptr<Search> SearchOP;

} // namespace motif_search

#endif // RNAMAKE_NEW_SEARCH_H
