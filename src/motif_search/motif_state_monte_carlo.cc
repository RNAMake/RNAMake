//
// Created by Joseph Yesselman on 11/4/17.
//

#include "base/log.hpp"
#include <fstream>

#include <motif_search/motif_state_monte_carlo.h>

namespace motif_search {

void MotifStateMonteCarlo::setup_options() {
  _options.add_option("sterics", false, base::OptionType::BOOL);
  _options.add_option("max_solutions", 1, base::OptionType::INT);
  _options.add_option("accept_score", 5, base::OptionType::FLOAT);
  _options.add_option("stages", 50, base::OptionType::INT);
  _options.add_option("steps", 500000, base::OptionType::INT);
  _options.lock_option_adding();

  update_var_options();
}

void MotifStateMonteCarlo::update_var_options() {
  _sterics = _options.get_bool("sterics");
  _max_solutions = _options.get_int("max_solutions");
  _accept_score = _options.get_float("accept_score");
  _stages = _options.get_int("stages");
  _steps = _options.get_int("steps");
}

void MotifStateMonteCarlo::setup(motif_data_structure::MotifStateGraphOP msg,
                                 int ni, int nj, int ei, int ej,
                                 bool target_an_aligned_end) {

  _msg = std::make_shared<motif_data_structure::MotifStateGraph>(*msg);
  _msg ->set_option_value("sterics", false);
  _ni = ni;
  _nj = nj;
  _ei = ei;
  _ej = ej;
  _target_an_aligned_end = target_an_aligned_end;

  _end = _msg->get_node(_nj)->data()->get_end_state(_ej);
  _end_flip =
      structure::BasepairStateOP(new structure::BasepairState(end_->copy()));
  _end_flip->flip();

  _org_num = _msg->size();

  auto i = -1;
  for (auto const &motif_states : _mses) {
    i++;
    auto ms = motif_states[_rng.randrange((int)motif_states.size() - 1)];
    if (i == 0) {
      _msg->add_state(ms, _ni, _ei);
    } else {
      _msg->add_state(ms);
    }
  }

  _stage = 0;
  _step = 0;
  _enumerating = false;
  _seen = StringIntMap();
}

MotifStateMonteCarloSolutionOP MotifStateMonteCarlo::next() {
  // do not call run() now, getting one solution at a time
  _enumerating = true;

  // initial search to not be clashing
  auto done = false;
  auto count = 0;
  while (!done) {
    count += 1;
    for (int i = 0; i < 10; i++) {
      perform_motif_swap_no_clash();
    }
    if (!_steric_clash(_msg)) {
      break;
    }

    if (count > 1000) {
      LOG_WARNING << "no suitable starting solutions!";
    }
  }

  auto score_num = 0;
  auto pos = 0;
  auto cur_score =
      get_score(_msg->last_node()->data()->cur_state->end_states()[1]);
  auto new_score = 0.0;
  auto best_score = cur_score;
  _mc.set_temperature(1.0);
  while (_stage < _stages) {
    while (_step < _steps) {
      _step += 1;
      new_score = perform_motif_swap(cur_score);
      // monte carlo move accepted
      if (new_score > 0) {
        cur_score = new_score;
      } else {
        continue;
      }

      if (new_score < best_score) {
        best_score = new_score;
      }

      // accept solution?
      if (new_score < _accept_score) {
        auto dist = 0.0;
        /*for(int k = org_num_; k < msg_->size(); k++) {
            for(auto const & b1 : msg_->get_node(k)->data()->cur_state->beads())
        { for(int l = 0; l < org_num_; l++) { for (auto const & b2 :
        msg_->get_node(l)->data()->cur_state->beads()) { dist = b1.distance(b2);
                        if(dist < 2.65) {
                            std::cout << "made it" << std::endl;
                        }
                    }
                }
            }
        }*/

        if (_seen_solution(_msg)) {
          continue;
        }
        auto mg = _msg->to_motif_graph();
        auto last_n = _msg->last_node()->index();
        auto nj_name = _msg->get_node(_nj)->data()->end_name(_ej);
        auto last_name = _msg->get_node(last_n)->data()->end_name(1);
        mg->add_connection(_nj, last_n, nj_name, last_name);
        return std::make_shared<MotifStateMonteCarloSolution>(mg, new_score);
      }
    }
    LOGI << "stage: " << _stage << " best_score: " << best_score
         << " cur_score: " << cur_score;

    int j = 0;
    for (int i = _org_num; i < _org_num + _mses.size(); i++) {
      bool found = 0;
      for (int k = 0; k < 100; k++) {
        auto new_ms = _mses[j][_rng.randrange((int)_mses[pos].size() - 1)];
        _msg->replace_state(i, new_ms);
        if (!_steric_clash(_msg)) {
          found = 1;
          break;
        }
      }
      if (!found) {
        LOG_WARNING << "could not find viable swap during heat up";
      }
      j++;
    }

    cur_score =
        get_score(_msg->last_node()->data()->cur_state->end_states()[1]);

    _step = 0;
    _stage += 1;
  }

  return MotifStateMonteCarloSolutionOP(nullptr);
}

MotifStateMonteCarloSolutionNewOP MotifStateMonteCarlo::next_state() {
  // do not call run() now, getting one solution at a time
  _enumerating = true;

  // initial search to not be clashing
  auto done = false;
  auto count = 0;
  while (!done) {
    count += 1;
    for (int i = 0; i < 10; i++) {
      perform_motif_swap_no_clash();
    }
    if (!_steric_clash(_msg)) {
      break;
    }

    if (count > 1000) {
      LOG_WARNING << "no suitable starting solutions!";
    }
  }

  auto score_num = 0;
  auto pos = 0;
  auto cur_score =
      get_score(_msg->last_node()->data()->cur_state->end_states()[1]);
  auto new_score = 0.0;
  auto best_score = cur_score;
  _mc.set_temperature(1.0);
  while (_stage < _stages) {
    while (_step < _steps) {
      _step += 1;
      new_score = perform_motif_swap(cur_score);
      // monte carlo move accepted
      if (new_score > 0) {
        cur_score = new_score;
      } else {
        continue;
      }

      if (new_score < best_score) {
        best_score = new_score;
      }

      // accept solution?
      if (new_score < _accept_score) {
        auto msg_copy =
            std::make_shared<motif_data_structure::MotifStateGraph>(*_msg);
        if (_seen_solution(_msg)) {
          continue;
        }
        auto last_n = _msg->last_node()->index();
        auto nj_name = _msg->get_node(_nj)->data()->end_name(_ej);
        auto last_name = _msg->get_node(last_n)->data()->end_name(1);
        msg_copy->add_connection(_nj, last_n, nj_name, last_name);
        return std::make_shared<MotifStateMonteCarloSolutionNew>(msg_copy,
                                                                 new_score);
      }
    }
    LOGI << "stage: " << _stage << " best_score: " << best_score
         << " cur_score: " << cur_score;

    int j = 0;
    for (int i = _org_num; i < _org_num + _mses.size(); i++) {
      bool found = 0;
      for (int k = 0; k < 100; k++) {
        auto new_ms = _mses[j][_rng.randrange((int)_mses[pos].size() - 1)];
        _msg->replace_state(i, new_ms);
        if (!_steric_clash(_msg)) {
          found = 1;
          break;
        }
      }
      if (!found) {
        LOG_WARNING << "could not find viable swap during heat up";
      }
      j++;
    }

    cur_score =
        get_score(_msg->last_node()->data()->cur_state->end_states()[1]);

    _step = 0;
    _stage += 1;
  }

  return MotifStateMonteCarloSolutionNewOP(nullptr);
}

void MotifStateMonteCarlo::start() { _enumerating = true; }

bool MotifStateMonteCarlo::finished() {
  if (!_enumerating) {
    throw std::runtime_error("finished() was called but was not enumerating?");
  }

  if (_seen.size() >= _max_solutions) {
    _enumerating = false;
    return true;
  } else {
    return false;
  }
}

void MotifStateMonteCarlo::run() {

  auto out = std::ofstream();
  out.open("default.scores");
  out << "design_num,design_score,motif_uses,topology\n";

  auto out_str = std::ofstream();
  out_str.open("default.out");

  auto seen = StringIntMap();

  auto score_num = 0;
  auto pos = 0;
  auto new_ms = motif::MotifStateOP(nullptr);
  auto last_ms = motif::MotifStateOP(nullptr);
  auto cur_score =
      get_score(_msg->last_node()->data()->cur_state->end_states()[1]);
  auto new_score = 0.0;
  auto best_score = cur_score;
  auto accept = 0;
  auto motif_used_string = String("");
  _mc.set_temperature(1.0);
  while (_stage < 100) {
    for (int i = 0; i < 500000; i++) {
      pos = _rng.randrange((int)_mses.size() - 1);
      new_ms = _mses[pos][_rng.randrange((int)_mses[pos].size() - 1)];
      last_ms = _msg->get_node(pos + _org_num)->data()->cur_state;
      _msg->replace_state(pos + _org_num, new_ms);
      new_score =
          get_score(_msg->last_node()->data()->cur_state->end_states()[1]);
      accept = _mc.accept(cur_score, new_score);
      if (best_score > cur_score) {
        best_score = cur_score;
      }
      if (accept) {
        if (_steric_clash(_msg)) {
          _msg->replace_state(pos + _org_num, last_ms);
          continue;
        }

        cur_score = new_score;
      } else {
        _msg->replace_state(pos + _org_num, last_ms);
        continue;
      }
      if (new_score < _accept_score) {
        int j = -1;
        motif_used_string = "";
        for (auto const &n : *_msg) {
          j++;
          if (j == 0) {
            continue;
          }
          motif_used_string += n->data()->name() + ";";
        }
        if (seen.find(motif_used_string) != seen.end()) {
          continue;
        }
        seen[motif_used_string] = 1;
        if (seen.size() % 10 == 0) {
          _msg->to_motif_graph()->to_pdb(
              "design." + std::to_string(seen.size() - 1) + ".pdb", 1);
          std::cout << "solutions: " << seen.size() << std::endl;
        }
        out << score_num << "," << cur_score << "," << motif_used_string << ","
            << std::endl;
        out_str << _msg->to_motif_graph()->to_str() << std::endl;
      }
    }

    // heat back up
    _mc.set_temperature(10.0);
    for (int i = 0; i < 1000; i++) {
      pos = _rng.randrange((int)_mses.size() - 1);
      new_ms = _mses[pos][_rng.randrange((int)_mses[pos].size() - 1)];
      last_ms = _msg->get_node(pos + _org_num)->data()->cur_state;
      _msg->replace_state(pos + _org_num, new_ms);
      new_score =
          get_score(_msg->last_node()->data()->cur_state->end_states()[1]);
      accept = _mc.accept(cur_score, new_score);
      if (accept) {
        if (_steric_clash(_msg)) {
          _msg->replace_state(pos + _org_num, last_ms);
          continue;
        }
        cur_score = new_score;
      } else {
        _msg->replace_state(pos + _org_num, last_ms);
      }
    }
    _mc.set_temperature(1.0);

    _stage += 1;
  }
}

double MotifStateMonteCarlo::get_score(structure::BasepairStateOP last_bp) {

  _score = last_bp->d().distance(_end->d());

  if (!_target_an_aligned_end) {
    _r_diff_flip = last_bp->r().difference(_end_flip->r());
  } else {
    _r_diff_flip = last_bp->r().difference(_end_flip->r());
  }
  _score += _r_diff_flip * 2;

  return _score;
}

bool MotifStateMonteCarlo::_steric_clash(
    motif_data_structure::MotifStateGraphOP msg) {

  // if (!sterics_) { return false; }
  auto clash = false;
  if (_using_lookup) {
    for (int i = _org_num; i < msg->size(); i++) {
      for (auto const &b : msg->get_node(i)->data()->cur_state->beads()) {
        clash = _lookup.clash(b);
        if (clash) {
          return true;
        }
      }
    }
  }
  return false;
}

float MotifStateMonteCarlo::perform_motif_swap(float cur_score) {
  auto new_ms = motif::MotifStateOP(nullptr);
  auto last_ms = motif::MotifStateOP(nullptr);
  auto pos = _rng.randrange((int)_mses.size() - 1);
  new_ms = _mses[pos][_rng.randrange((int)_mses[pos].size() - 1)];
  last_ms = _msg->get_node(pos + _org_num)->data()->cur_state;
  _msg->replace_state(pos + _org_num, new_ms);
  auto new_score =
      get_score(_msg->last_node()->data()->cur_state->end_states()[1]);

  // should we accept?
  auto accept = _mc.accept(cur_score, new_score);

  if (accept) {
    if (_steric_clash(_msg)) {
      _msg->replace_state(pos + _org_num, last_ms);
      return -1;
    }
    return new_score;
  } else {
    _msg->replace_state(pos + _org_num, last_ms);
    return -1;
  }
}

void MotifStateMonteCarlo::perform_motif_swap_no_clash() {
  auto new_ms = motif::MotifStateOP(nullptr);
  auto last_ms = motif::MotifStateOP(nullptr);
  auto pos = _rng.randrange((int)_mses.size() - 1);
  new_ms = _mses[pos][_rng.randrange((int)_mses[pos].size() - 1)];
  _msg->replace_state(pos + _org_num, new_ms);
}

bool MotifStateMonteCarlo::_seen_solution(
    motif_data_structure::MotifStateGraphOP msg) {

  int j = -1;
  auto motif_used_string = String("");
  for (auto const &n : *_msg) {
    j++;
    if (j == 0) {
      continue;
    }
    if (n->data()->name().substr(0, 10) == "HELIX.FLEX") {
      auto spl = base::string::split(n->data()->name(), ".");
      motif_used_string += spl[0] + "." + spl[1] + "." + spl[2] + ";";
    } else {
      motif_used_string += n->data()->name() + ";";
    }
  }
  if (_seen.find(motif_used_string) != _seen.end()) {
    return true;
  } else {
    _seen[motif_used_string] = 1;
    return false;
  }
}

} // namespace motif_search
