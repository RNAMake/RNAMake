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

  /// @brief - constructors
  inline explicit Residue(util::Beads& beads): _beads(std::move(beads)) {}
  // constructor
  Residue(const Residue &) = default;
  // destructor
  ~Residue() = default;

public: // getters
  typedef util::Beads::const_iterator const_iterator;

  /// @brief - gets the beginning (?) of the
  [[nodiscard]] const_iterator begin() const noexcept { return _beads.begin(); }

  /// @brief -
  [[nodiscard]] const_iterator end() const noexcept { return _beads.end(); }

  /// @brief -
  [[nodiscard]] inline size_t num_of_beads() const { return _beads.size(); }

public:
  /// @brief - moves the residue by the appropriate vector
  void move(const math::Vector3 &p) {}


private:
  util::Beads _beads;
};

typedef std::vector<Residue> Residues;

} // namespace structure::state

#endif // RNAMAKE_SRC_STRUCTURE_STATE_RESIDUE_HPP_
