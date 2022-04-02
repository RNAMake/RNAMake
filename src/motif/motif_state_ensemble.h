//
//  motif_state_ensemble.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_ensemble__
#define __RNAMake__motif_state_ensemble__

#include <stdio.h>

#include <algorithm>

#include "motif/motif_state.h"

namespace motif {

/*
 * Exception for motif state ensemble
 */
class MotifStateEnsembleException : public std::runtime_error {
 public:
  /**
   * Standard constructor for MotifStateEnsembleException
   * @param   message   Error message for motif state ensembles
   */
  MotifStateEnsembleException(String const &message)
      : std::runtime_error(message) {}
};

struct MotifStateEnsembleMember {
  inline MotifStateEnsembleMember(MotifStateOP const &nmotif_state,
                                  float const &nenergy)
      : motif_state(nmotif_state), energy(nenergy) {}

  inline MotifStateEnsembleMember(MotifStateEnsembleMember const &msem)
      : motif_state(std::make_shared<MotifState>(*msem.motif_state)),
        energy(msem.energy) {}

  inline String to_str() {
    return motif_state->to_str() + "#" + std::to_string(energy);
  }

  MotifStateOP motif_state;
  float energy;
};

typedef std::shared_ptr<MotifStateEnsembleMember> MotifStateEnsembleMemberOP;
typedef std::vector<MotifStateEnsembleMemberOP> MotifStateEnsembleMemberOPs;

// to allow sorting of members
struct MotifStateEnsembleMember_LessThanKey {
  inline bool operator()(MotifStateEnsembleMemberOP const &mem1,
                         MotifStateEnsembleMemberOP const &mem2) {
    return mem1->energy < mem2->energy;
  }
};

class MotifStateEnsemble {
 public:
  MotifStateEnsemble(MotifStateOP const &ms)
      : id_(ms->end_ids()[0]),
        block_end_add_(ms->block_end_add()),
        members_(MotifStateEnsembleMemberOPs()) {
    members_.push_back(std::make_shared<MotifStateEnsembleMember>(ms, 1));
  }

  MotifStateEnsemble(MotifStateEnsemble const &mse)
      : id_(mse.id_),
        block_end_add_(mse.block_end_add_),
        members_(MotifStateEnsembleMemberOPs()) {
    for (auto const &mem : mse.members_) {
      auto new_mem = std::make_shared<MotifStateEnsembleMember>(*mem);
      members_.push_back(new_mem);
    }
  }

  MotifStateEnsemble(MotifStateOPs const &, Floats const &);

  MotifStateEnsemble(String const &);

  ~MotifStateEnsemble() {}

 public:
  typedef typename MotifStateEnsembleMemberOPs::const_iterator const_iterator;
  typedef typename MotifStateEnsembleMemberOPs::iterator iterator;

  iterator begin() { return members_.begin(); }
  iterator end() { return members_.end(); }

  const_iterator begin() const noexcept { return members_.begin(); }
  const_iterator end() const noexcept { return members_.end(); }

 public:
  String to_str();

  inline int num_end_states() const {
    return (int)members_[0]->motif_state->end_names().size();
  }

  inline MotifStateOP const &most_populated() const {
    return members_[0]->motif_state;
  }

  inline MotifStateEnsembleMemberOP const &get_member(int i) const {
    if (i >= members_.size()) {
      throw MotifStateEnsembleException("cannot get member " +
                                        std::to_string(i) +
                                        " it does not exist, size of members is"
                                        " " +
                                        std::to_string(members_.size()));
    }

    return members_[i];
  }

  inline int member_index(MotifStateEnsembleMemberOP const &mem) {
    return std::find(members_.begin(), members_.end(), mem) - members_.begin();
  }

  inline MotifStateEnsembleMemberOPs const &members() { return members_; }

 public:
  size_t size() const { return members_.size(); }

 public:
  inline int block_end_add() { return block_end_add_; }

  inline String id() { return id_; }

 private:
  String id_;
  int block_end_add_;
  MotifStateEnsembleMemberOPs members_;
};

typedef std::shared_ptr<MotifStateEnsemble> MotifStateEnsembleOP;
typedef std::vector<MotifStateEnsembleOP> MotifStateEnsembleOPs;

}  // namespace motif

#endif /* defined(__RNAMake__motif_state_ensemble__) */
