//
// Created by Joseph Yesselman on 2019-03-30.
//

#include "motif_search/exhaustive/motif_state_enumerator.h"

namespace motif_search {
namespace exhaustive {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// constructors
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

MotifStateEnumerator::MotifStateEnumerator(
    motif_search::SolutionToplogy sol_toplogy) {
  indices_ = Indexes(sol_toplogy.size());
  maxes_ = Indexes(sol_toplogy.size());
  motif_states_ = std::vector<motif::MotifStateOPs>(sol_toplogy.size());
  current_ = motif::MotifStateOPs(sol_toplogy.size());
  end_ = 0;
  size_limit_ = 99999;
  int i = 0;
  auto m_states = motif::MotifStateOPs();
  for (auto const &n : sol_toplogy) {
    m_states.resize(0);
    for (auto const &mem : *n->data()) {
      if (mem->motif_state->name().substr(0, 10) == "HELIX.FLEX") {
        if (mem->motif_state->size() > 14) {
          continue;
        }
      }
      mem->motif_state->new_uuids();
      m_states.push_back(
          std::make_shared<motif::MotifState>(*mem->motif_state));
    }
    motif_states_[i] = m_states;
    maxes_[i] = m_states.size();
    indices_[i] = 0;
    i++;
  }
  if (motif_states_.size() == 0) {
    throw std::runtime_error("must supply a topology with atleast one node");
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main interface
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MotifStateEnumerator::start(structure::BasepairStateOP start_bp) {
  start_bp_ = start_bp;
  update_ = 0;
  end_ = 0;

  for (auto &index : indices_) {
    index = 0;
  }
}

bool MotifStateEnumerator::finished() { return end_; }

void MotifStateEnumerator::next() {
  auto done = false;
  while (!done) {
    _iterate();
    if (_within_size_limit()) {
      break;
    }
  }
}

motif::MotifStateOP MotifStateEnumerator::top_state() {
  if (!end_ || !updated_) {
    _update_current_states();
  }
  return motif_states_.back()[indices_.back()];
}

motif::MotifStateOPs const &MotifStateEnumerator::all_states() {
  if (!end_ || !updated_) {
    _update_current_states();
  }
  int j = 0;
  for (auto const &v : motif_states_) {
    current_[j] = v[indices_[j]];
    j++;
  }
  return current_;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MotifStateEnumerator::_update_current_states() {
  updated_ = true;
  if (update_ == 0) {
    aligner_.get_aligned_motif_state(start_bp_, motif_states_[0][indices_[0]]);
    update_ += 1;
  }
  for (int j = update_; j < indices_.size(); j++) {
    aligner_.get_aligned_motif_state(
        motif_states_[j - 1][indices_[j - 1]]->end_states()[1],
        motif_states_[j][indices_[j]]);
  }
}

void MotifStateEnumerator::_iterate() {
  bool not_updated = true;
  if (updated_ == false) {
    not_updated = false;
  }
  updated_ = false;
  int i = (int)indices_.size() - 1;
  while (i > -1) {
    if (not_updated) {
      update_ = i;
    } else {
      if (updated_ > i) {
        updated_ = i;
      }
    }
    indices_[i]++;
    if (indices_[i] == maxes_[i]) {
      if (i == 0) {
        end_ = 1;
        indices_[i]--;
        break;
      }
      indices_[i] = 0;
      i--;
      continue;
    }
    break;
  }
}

bool MotifStateEnumerator::_within_size_limit() {
  auto total_size = -2;
  auto j = 0;
  for (auto const &v : motif_states_) {
    total_size += v[indices_[j]]->size() - 2;
    j += 1;
  }
  if (total_size > size_limit_) {
    return false;
  }
  return true;
}

} // namespace exhaustive
} // namespace motif_search
