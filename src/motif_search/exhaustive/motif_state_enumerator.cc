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
  _indices = Indexes(sol_toplogy.size());
  _maxes = Indexes(sol_toplogy.size());
  _motif_states = std::vector<motif::MotifStateOPs>(sol_toplogy.size());
  _current = motif::MotifStateOPs(sol_toplogy.size());
  _end = 0;
  _size_limit = 99999;
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
    _motif_states[i] = m_states;
    _maxes[i] = m_states.size();
    _indices[i] = 0;
    i++;
  }
  if (_motif_states.size() == 0) {
    throw std::runtime_error("must supply a topology with atleast one node");
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main interface
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MotifStateEnumerator::start(structure::BasepairStateOP start_bp) {
  _start_bp = start_bp;
  _update = 0;
  _end = 0;

  for (auto &index : _indices) {
    index = 0;
  }
}

bool MotifStateEnumerator::finished() { return _end; }

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
  if (!_end || !_updated) {
    _update_current_states();
  }
  return _motif_states.back()[_indices.back()];
}

motif::MotifStateOPs const &MotifStateEnumerator::all_states() {
  if (!_end || !_updated) {
    _update_current_states();
  }
  int j = 0;
  for (auto const &v : _motif_states) {
    _current[j] = v[_indices[j]];
    j++;
  }
  return _current;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MotifStateEnumerator::_update_current_states() {
  _updated = true;
  if (_update == 0) {
    _aligner.get_aligned_motif_state(_start_bp, _motif_states[0][_indices[0]]);
    _update += 1;
  }
  for (int j = _update; j < _indices.size(); j++) {
    _aligner.get_aligned_motif_state(
        _motif_states[j - 1][_indices[j - 1]]->end_states()[1],
        _motif_states[j][_indices[j]]);
  }
}

void MotifStateEnumerator::_iterate() {
  bool not_updated = true;
  if (_updated == false) {
    not_updated = false;
  }
  _updated = false;
  int i = (int)_indices.size() - 1;
  while (i > -1) {
    if (not_updated) {
      _update = i;
    } else {
      if (_updated > i) {
        _updated = i;
      }
    }
    _indices[i]++;
    if (_indices[i] == _maxes[i]) {
      if (i == 0) {
        _end = 1;
        _indices[i]--;
        break;
      }
      _indices[i] = 0;
      i--;
      continue;
    }
    break;
  }
}

bool MotifStateEnumerator::_within_size_limit() {
  auto total_size = -2;
  auto j = 0;
  for (auto const &v : _motif_states) {
    total_size += v[_indices[j]]->size() - 2;
    j += 1;
  }
  if (total_size > _size_limit) {
    return false;
  }
  return true;
}

} // namespace exhaustive
} // namespace motif_search
