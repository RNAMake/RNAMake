//
// Created by Joseph Yesselman on 3/23/19.
//

#ifndef RNAMAKE_NEW_MONTE_CARLO_SEARCH_H
#define RNAMAKE_NEW_MONTE_CARLO_SEARCH_H

#include <base/types.hpp>
//#include <base/option.h>
//#include <motif/motif_state_aligner.h>
#include <motif_data_structure/motif_state_graph.hpp>
#include <motif_search/monte_carlo/scorer.h>
#include <motif_search/search.h>
#include <motif_search/solution_filter.h>
#include <motif_search/solution_topology.h>
#include <util/monte_carlo.h>

namespace motif_search {
namespace monte_carlo {

class Move {
public:
  Move(String const &name) : _name(name), _mc(util::MonteCarlo(3.0)) {}

  virtual ~Move() {}

public:
  virtual bool apply(motif_data_structure::MotifStateGraphOP, float) = 0;

  virtual float score() = 0;

  virtual void undo(motif_data_structure::MotifStateGraphOP) = 0;

public:
  void set_temperature(float temp) { _mc.set_temperature(temp); }

  void scale_temperature(float scale) { _mc.scale_temperature(scale); }

protected:
  String _name;
  util::MonteCarlo _mc;
};

typedef std::shared_ptr<Move> MoveOP;
typedef std::vector<MoveOP> MoveOPs;

class MotifSwapMove : public Move {
public:
  MotifSwapMove(ScorerOP scorer, SolutionToplogy const &sol_top)
      : Move("MotifSwap"), _scorer(scorer), _sol_top(sol_top) {
    this->_mc = util::MonteCarlo(1.0f);
    _rng = util::RandomNumberGenerator();
  }

public:
  bool apply(motif_data_structure::MotifStateGraphOP msg, float current_score) {
    _pos = _rng.randrange((int)_sol_top.size() - 1);
    _new_ms = _sol_top.get_motif_state(_pos);
    _last_ms = msg->get_node(_pos + 1)->data()->cur_state;
    msg->replace_state(_pos + 1, _new_ms);
    _new_score =
        _scorer->score(*msg->last_node()->data()->cur_state->end_states()[1]);
    _accept = _mc.accept(current_score, _new_score);
    if (_accept) {
      return true;
    } else {
      undo(msg);
      return false;
    }

    return true;
  }

  float score() { return _new_score; }

  void undo(motif_data_structure::MotifStateGraphOP msg) {
    msg->replace_state(_pos + 1, _last_ms);
  }

private:
  ScorerOP _scorer;
  SolutionToplogy _sol_top;
  float _new_score;
  int _accept, _pos;
  util::RandomNumberGenerator _rng;
  motif::MotifStateOP _new_ms, _last_ms;
};

class MoveSet {};

class Search : public motif_search::Search {
public:
  struct Parameters {
    float accept_score = 10.0f;
    int max_size = 100000;
  };

public:
  Search(ScorerOP scorer, SolutionToplogy const &sol_top,
         SolutionFilterOP filter)
      : motif_search::Search("monte_carlo"), _scorer(scorer->clone()),
        _sol_top(sol_top), _filter(filter->clone()) {
    _finished = false;
    setup_options();
    update_var_options();
  }

  ~Search() {}

  motif_search::Search *clone() const { return new Search(*this); };

public:
  virtual void setup(ProblemOP p) {
    _msg = _sol_top.initialize_solution(p->start);
    _scorer->set_target(p->end, p->target_an_aligned_end);
    _lookup = p->lookup;
    if (_lookup != nullptr) {
      _using_lookup = true;
    }
    _stages = 50;
    _steps = 500000;
    _stage = 0;
  }

  virtual void
  // TODO implement or change interface
  start() {}

  virtual bool
  // TODO implement or change interface
  finished() {
    return false;
  }

  virtual SolutionOP next() {
    auto cur_score =
        _scorer->score(*_msg->last_node()->data()->cur_state->end_states()[1]);
    auto new_score = 0.0;
    auto best_score = cur_score;
    auto mover = std::make_shared<MotifSwapMove>(_scorer, _sol_top);
    auto hot_mover = std::make_shared<MotifSwapMove>(_scorer, _sol_top);
    auto min_mover = std::make_shared<MotifSwapMove>(_scorer, _sol_top);
    hot_mover->set_temperature(100.0f);
    min_mover->set_temperature(0.1f);
    auto temp = 4.5;
    mover->set_temperature(temp);

    // TODO add temperature adjustment to get to 0.235 acceptance
    // TODO add mnimization step
    // TODO add sterics
    auto accept = false;
    auto accepted_steps = 0.0;
    while (_stage < _stages) {
      accepted_steps = 0.0;
      auto round_best_score = 100000.0f;
      auto best_sol = motif_data_structure::MotifStateGraphOP(nullptr);
      while (_step < _steps) {
        _step += 1;
        accept = mover->apply(_msg, cur_score);
        if (!accept) {
          continue;
        }
        accepted_steps += 1;
        cur_score = mover->score();
        if (cur_score < _parameters.accept_score) {
          _get_solution_motif_names(_msg);
          if (!_filter->accept(_motif_names)) {
            continue;
          }
          auto sol_msg = _get_solution_msg();
          LOG_DEBUG << "found a solution: " << cur_score;
          return std::make_shared<Solution>(sol_msg, cur_score);
        }
        if (cur_score < best_score) {
          best_score = cur_score;
        }
        if (cur_score < round_best_score) {
          round_best_score = cur_score;
        }
      }
      auto accept_ratio = (float)(accepted_steps / _steps);
      LOG_DEBUG << "stage: " << _stage << " best_score: " << best_score
                << " acceptance: " << accept_ratio << " temp: " << temp;
      LOG_DEBUG << "round best score: " << round_best_score;
      // heatup
      for (int i = 0; i < 1000; i++) {
        accept = hot_mover->apply(_msg, cur_score);
        if (accept) {
          cur_score = hot_mover->score();
        }
      }
      // reset mover
      mover = std::make_shared<MotifSwapMove>(_scorer, _sol_top);
      auto diff = (accept_ratio - 0.235) * 10;
      // temp = temp - diff;
      // mover->set_temperature(temp);

      _step = 0;
      _stage += 1;
    }

    LOG_DEBUG << "exiting monte carlo search";

    return SolutionOP(nullptr);
  }

private:
  bool _steric_clash(motif_data_structure::MotifStateGraphOP msg) {
    auto clash = false;
    if (_using_lookup) {
      for (auto const &n : *msg) {
        for (auto const &b : n->data()->cur_state->beads()) {
          clash = _lookup->clash(b);
          if (clash) {
            return true;
          }
        }
      }
    }
    return false;
  }

  motif_data_structure::MotifStateGraphOP _get_solution_msg() {
    auto new_msg = std::make_shared<motif_data_structure::MotifStateGraph>();
    new_msg->set_option_value("sterics", false);
    for (auto const &n : *_msg) {
      if (n->index() == 0) {
        continue;
      }
      if (n->index() == 1) {
        new_msg->add_state(
            std::make_shared<motif::MotifState>(*n->data()->cur_state));
      } else {
        new_msg->add_state(
            std::make_shared<motif::MotifState>(*n->data()->cur_state), -1,
            n->parent_end_index());
      }
    }
    return new_msg;
  }

  void _get_solution_motif_names(motif_data_structure::MotifStateGraphOP msg) {
    _motif_names.resize(0);
    for (auto const &n : *msg) {
      _motif_names.push_back(n->data()->name());
    }
  }

protected:
  void setup_options() {
    _options.add_option("sterics", true, base::OptionType::BOOL);
    _options.add_option("min_size", 0, base::OptionType::INT);
    _options.add_option("max_size", 1000000, base::OptionType::INT);
    _options.add_option("max_solutions", 1, base::OptionType::INT);
    _options.add_option("accept_score", 10.0f, base::OptionType::FLOAT);
    _options.add_option("return_best", false, base::OptionType::BOOL);
    _options.lock_option_adding();
  }

  void update_var_options() {
    _parameters.accept_score = _options.get_float("accept_score");
    _parameters.max_size = _options.get_int("max_size");
  }

private:
  Parameters _parameters;
  ScorerOP _scorer;
  SolutionToplogy _sol_top;
  SolutionFilterOP _filter;
  util::StericLookupNewOP _lookup;
  motif_data_structure::MotifStateGraphOP _msg;
  float _score;
  int _stages, _stage;
  int _steps, _step;
  bool _finished, _using_lookup;
  Strings _motif_names;
};

} // namespace monte_carlo
} // namespace motif_search

#endif // RNAMAKE_NEW_MONTE_CARLO_SEARCH_H
