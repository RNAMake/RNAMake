//
// Created by Joe Yesselman on 6/28/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_STATE_RESIDUE_HPP_
#define RNAMAKE_SRC_STRUCTURE_STATE_RESIDUE_HPP_

#include <math/rotation.hpp>
#include <util/bead.h>

namespace structure::state {

class Residue {
public:
  inline explicit Residue(util::Beads& beads): _beads(std::move(beads)) {}

  Residue(const Residue &) = default;
  
  ~Residue() = default;
public:
  typedef util::Beads::const_iterator const_iterator;

  [[nodiscard]] const_iterator begin() const noexcept { return _beads.begin(); }

  [[nodiscard]] const_iterator end() const noexcept { return _beads.end(); }

  [[nodiscard]] inline size_t num_of_beads() const {return _beads.size(); }

public:
  void move(const math::Vector3 &p) {}

  void transform(const math::RotandTrans &rt) {}

private:
  util::Beads _beads;
};

typedef std::vector<Residue> Residues;

} // namespace structure::state

#endif // RNAMAKE_SRC_STRUCTURE_STATE_RESIDUE_HPP_
