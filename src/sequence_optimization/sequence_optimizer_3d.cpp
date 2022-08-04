//
//  sequence_optimizer_3d.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/30/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "sequence_optimization/sequence_optimizer_3d.hpp"
#include "secondary_structure/util.h"
#include "util/monte_carlo.h"
#include <base/log.hpp>

/*
namespace sequence_optimization {

SequenceOptimizer3D::SequenceOptimizer3D()
    : _eterna_scorer(eternabot::Scorer()), _scorer(nullptr),
      _rng(util::RandomNumberGenerator()) {
  _possible_bps =
      std::vector<Strings>({{"A", "U"}, {"U", "A"}, {"G", "C"}, {"C", "G"}});

  _disallowed_sequences = Strings();
  _disallowed_sequences.push_back(String("AAAA"));
  _disallowed_sequences.push_back(String("CCCC"));
  _disallowed_sequences.push_back(String("GGGG"));
  _disallowed_sequences.push_back(String("UUUU"));

  _disallowed_res_types_sequences =
      std::vector<secondary_structure::ResTypes>();
  for (auto const &seq : _disallowed_sequences) {
    auto disallowed_types = secondary_structure::ResTypes();
    secondary_structure::get_res_types_from_sequence(seq, disallowed_types);
    _disallowed_res_types_sequences.push_back(disallowed_types);
  }
  _current_violations = Indexes(_disallowed_res_types_sequences.size());
  _next_violations = Indexes(_disallowed_res_types_sequences.size());

  setup_options();
}

// get_optimized_sequence helpers
// //////////////////////////////////////////////////////////////////

void SequenceOptimizer3D::_update_designable_bp(
    DesignableBPOP const &d_bp, motif_data_structure::MotifStateGraphOP &msg,
    secondary_structure::PoseOP &ss) {

  if (d_bp->m_id_bot != nullptr) {
    ss->update_motif(*d_bp->m_id_bot);
    auto m = ss->motif(*d_bp->m_id_bot);
    auto n = msg->get_node(*d_bp->m_id_bot);
    msg->replace_state(n->index(), resources::Manager::instance().bp_step_state(
                                       m->end_ids()[0]));
  }
  if (d_bp->m_id_top != nullptr) {
    ss->update_motif(*d_bp->m_id_top);
    auto m = ss->motif(*d_bp->m_id_top);
    auto n = msg->get_node(*d_bp->m_id_top);
    msg->replace_state(n->index(), resources::Manager::instance().bp_step_state(
                                       m->end_ids()[0]));
  }
}

String SequenceOptimizer3D::_validate_sequence(
    motif_data_structure::MotifStateGraphOP const &msg,
    secondary_structure::PoseOP const &ss) {

  auto s1 = msg->to_motif_graph()->secondary_structure()->sequence();
  auto s2 = ss->sequence();

  for (int j = 0; j < s2.length(); j++) {
    if (s1[j] != s2[j]) {
      // std::cout << s1[j] << " " << s2[j] << std::endl;
      throw std::runtime_error(
          "sequences are out of sync: something went really wrong in sequence "
          "optimization");
    }
  }

  return s1;
}

void SequenceOptimizer3D::find_seq_violations(secondary_structure::PoseOP ss,
                                              Indexes &violations) {
  auto pos = 0;
  for (int i = 0; i < violations.size(); i++) {
    violations[i] = 0;
  }

  // TODO go back and see if this works!!, maybe right a sub class to test
  // seperately
  for (auto const &c : ss->chains()) {
    auto &res = c->residues();
    for (int i = 0; i < res.size(); i++) {
      for (int j = 0; j < _disallowed_res_types_sequences.size(); j++) {
        if (i < _disallowed_res_types_sequences[j].size()) {
          continue;
        }
        auto match = true;
        pos = i - _disallowed_res_types_sequences[j].size();
        for (auto const &e : _disallowed_res_types_sequences[j]) {
          if (res[pos]->res_type() != e) {
            match = false;
            break;
          }
          pos += 1;
        }
        if (match) {
          violations[j] += 1;
        }
      }
    }
  }
}

int SequenceOptimizer3D::find_gc_helix_stretches(
    secondary_structure::PoseOP ss) {
  int count = 0;
  int violations = 0;
  for (auto const &h : ss->helices()) {
    count = 0;
    for (auto const &r : h->chains()[0]->residues()) {
      if (r->res_type() == secondary_structure::ResType::GUA ||
          r->res_type() == secondary_structure::ResType::CYT) {
        count += 1;
      } else {
        count = 0;
      }
      if (count > 3) {
        violations += 1;
        break;
      }
    }
  }
  return violations;
}

SequenceOptimizer3D::DesignableBPOPs
SequenceOptimizer3D::_get_designable_bps(secondary_structure::PoseOP &ss) {

  auto designable_bps = DesignableBPOPs();
  for (auto const &bp : ss->basepairs()) {
    auto bp_name = bp->res1()->name() + bp->res2()->name();
    if (bp_name == "NN") {
      designable_bps.push_back(std::make_shared<DesignableBP>(bp));
    }
  }

  find_seq_violations(ss, _current_violations);
  _current_gc_stretches = find_gc_helix_stretches(ss);
  int count = 0;

  for (auto const &d_bp : designable_bps) {
    for (auto const &m : ss->motifs()) {
      if (m->name() != "HELIX.IDEAL") {
        continue;
      }
      if (m->ends()[0] == d_bp->bp) {
        d_bp->m_id_bot = std::make_shared<util::Uuid>(m->id());
      }
      if (m->ends()[1] == d_bp->bp) {
        d_bp->m_id_top = std::make_shared<util::Uuid>(m->id());
      }
    }

    auto state = _possible_bps[_rng.randrange(_possible_bps.size())];
    d_bp->bp->res1()->name(state[0]);
    d_bp->bp->res2()->name(state[1]);
    find_seq_violations(ss, _next_violations);
    _next_gc_stretches = find_gc_helix_stretches(ss);
    while (new_seq_violations()) {
      auto state = _possible_bps[_rng.randrange(_possible_bps.size())];
      d_bp->bp->res1()->name(state[0]);
      d_bp->bp->res2()->name(state[1]);
      find_seq_violations(ss, _next_violations);
      _next_gc_stretches = find_gc_helix_stretches(ss);
      count++;
      if (count > 100) {
        _current_violations = _next_violations;
        _current_gc_stretches = _next_gc_stretches;
      }
    }
  }
  for (auto &m : ss->motifs()) {
    ss->update_motif(m->id());
  }

  return designable_bps;
}

void SequenceOptimizer3D::_initiate_sequence_in_msg(
    motif_data_structure::MotifStateGraphOP &msg,
    secondary_structure::PoseOP const &ss) {

  for (auto const &m : ss->motifs()) {
    if (m->name() != "HELIX.IDEAL") {
      continue;
    }

    auto name = String("");
    auto n_msg = msg->get_node(m->id());

    msg->replace_state(
        n_msg->index(),
        resources::Manager::instance().bp_step_state(m->end_ids()[0]));
  }
}

// optimizing sequence
// /////////////////////////////////////////////////////////////////////////////

OptimizedSequenceOPs SequenceOptimizer3D::get_optimized_sequences(
    motif_data_structure::MotifGraphOP const &mg) {
  update_var_options();

  if (_scorer == nullptr) {
    throw std::runtime_error(
        "cannot run get_optimized_sequences without scorer, either supply here "
        "or use set_scorer");
  }

  auto sols = OptimizedSequenceOPs();
  auto ss = mg->designable_secondary_structure();
  auto msg = std::make_shared<motif_data_structure::MotifStateGraph>(mg);

  //if(ss->chains().size() > 1) {
  //    throw std::runtime_error(
  //        "cannot perform sequence optmization with more than one chain");
  //}

  auto designable_bps = _get_designable_bps(ss);
  _initiate_sequence_in_msg(msg, ss);
  _eterna_scorer.setup(ss);

  auto last_score = _scorer->score(msg);
  auto new_score = 0.0f, eterna_score = 0.0f;

  // std::cout << designable_bps.size() << std::endl;
  auto mc = util::MonteCarlo(1.0f);

  auto best = 1000.0f;
  auto best_seq = String("");

  int i = -1;
  auto d_bp = DesignableBPOP(nullptr);
  auto new_bp_state = Strings();
  while (i < _steps) {
    i++;
    d_bp = designable_bps[_rng.randrange(designable_bps.size())];
    new_bp_state = _possible_bps[_rng.randrange(_possible_bps.size())];
    d_bp->update_state(new_bp_state);

    find_seq_violations(ss, _next_violations);
    _next_gc_stretches = find_gc_helix_stretches(ss);

    if (new_seq_violations()) {
      d_bp->revert_state();
      continue;
    }

    _update_designable_bp(d_bp, msg, ss);

    new_score = _scorer->score(msg);

    if (mc.accept(last_score, new_score)) {
      last_score = new_score;
    } else {
      d_bp->revert_state();
      _update_designable_bp(d_bp, msg, ss);
      continue;
    }

    if (best > new_score) {
      best = new_score;
      best_seq = ss->sequence();
      LOG_VERBOSE << "best_score=" << best;
    }

    if (_cutoff < new_score) {
      continue;
    }

    eterna_score = _eterna_scorer.score_secondary_structure(ss);
    if (eterna_score > _eterna_cutoff) {
      auto seq = _validate_sequence(msg, ss);
      sols.push_back(std::make_shared<OptimizedSequence>(
          OptimizedSequence{seq, new_score, eterna_score}));
      LOG_DEBUG << "found solution! score=" << new_score;
      LOG_DEBUG << seq;
      if (sols.size() >= _solutions) {
        return sols;
      }
    }
  }

  if (sols.size() == 0 && _return_lowest) {
    LOG_DEBUG << "could not reach cutoff but returning best solution=" << best;
    LOG_DEBUG << best_seq;
    ss->replace_sequence(best_seq);
    eterna_score = _eterna_scorer.score_secondary_structure(ss);
    sols.push_back(std::make_shared<OptimizedSequence>(
        OptimizedSequence{best_seq, best, eterna_score}));
  }

  return sols;
}

motif_data_structure::MotifGraphOP SequenceOptimizer3D::get_optimized_mg(
    motif_data_structure::MotifGraphOP const &mg) {
  update_var_options();

  if (_scorer == nullptr) {
    throw std::runtime_error(
        "cannot run get_optimized_sequences without scorer, either supply here "
        "or use set_scorer");
  }

  auto ss = mg->designable_secondary_structure();
  auto msg = std::make_shared<motif_data_structure::MotifStateGraph>(mg);

  if (ss->chains().size() > 1) {
    throw std::runtime_error(
        "cannot perform sequence optmization with more than one chain");
  }

  auto designable_bps = _get_designable_bps(ss);
  _initiate_sequence_in_msg(msg, ss);
  _eterna_scorer.setup(ss);

  auto last_score = _scorer->score(msg);
  auto new_score = 0.0f, eterna_score = 0.0f;

  auto mc = util::MonteCarlo(1.0f);
  auto best = 1000.0f;
  auto best_seq = String("");

  int i = -1;
  auto d_bp = DesignableBPOP(nullptr);
  auto new_bp_state = Strings();
  auto best_msg = std::make_shared<motif_data_structure::MotifStateGraph>();
  while (i < _steps) {
    i++;

    d_bp = designable_bps[_rng.randrange(designable_bps.size())];
    new_bp_state = _possible_bps[_rng.randrange(_possible_bps.size())];
    d_bp->update_state(new_bp_state);

    _update_designable_bp(d_bp, msg, ss);

    new_score = _scorer->score(msg);

    if (mc.accept(last_score, new_score)) {
      last_score = new_score;
    } else {
      d_bp->revert_state();
      _update_designable_bp(d_bp, msg, ss);
      continue;
    }

    if (best > new_score) {
      best = new_score;
      best_msg = std::make_shared<motif_data_structure::MotifStateGraph>(*msg);
      if (_verbose) {
        std::cout << "SEQUENCE OPTIMIZER: best_score=" << best << std::endl;
      }
    }

    if (_cutoff < new_score) {
      continue;
    }

    eterna_score = _eterna_scorer.score_secondary_structure(ss);
    if (eterna_score > _eterna_cutoff) {

      auto seq = _validate_sequence(msg, ss);
      return msg->to_motif_graph();

      if (_verbose) {
        std::cout << "SEQUENCE OPTIMIZER: found solution! score=" << new_score;
        std::cout << " eterna_score=" << eterna_score << std::endl;
      }
    }
  }

  return best_msg->to_motif_graph();
}

// options
// //////////////////////////////////////////////////////////////////////////////////////////

void SequenceOptimizer3D::setup_options() {
  _options.add_option("cutoff", 5.0f, base::OptionType::FLOAT);
  _options.add_option("solutions", 1, base::OptionType::INT);
  _options.add_option("eterna_cutoff", -1.0f, base::OptionType::FLOAT);
  _options.add_option("verbose", false, base::OptionType::BOOL);
  _options.add_option("return_lowest", true, base::OptionType::BOOL);
  _options.add_option("steps", 10000, base::OptionType::INT);
  _options.lock_option_adding();
  update_var_options();
}

void SequenceOptimizer3D::update_var_options() {
  _cutoff = _options.get_float("cutoff");
  _solutions = _options.get_int("solutions");
  _eterna_cutoff = _options.get_float("eterna_cutoff");
  _verbose = _options.get_bool("verbose");
  _return_lowest = _options.get_bool("return_lowest");
  _steps = _options.get_int("steps");
}

} // namespace sequence_optimization
*/
