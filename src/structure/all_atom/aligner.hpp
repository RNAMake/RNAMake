//
// Created by Joe Yesselman on 7/10/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_ALL_ATOM_ALIGNER_HPP_
#define RNAMAKE_SRC_STRUCTURE_ALL_ATOM_ALIGNER_HPP_

#include <structure/all_atom/segment.hpp>

namespace structure::all_atom {
class Aligner {
public:

  /// @brief - places a segment inside the chain at a certain index
  void align(const Segment &, Segment &, Index);

  // Segment get_aligned();

private:
  math::Matrix3x3 _rot;
};

} // namespace structure::all_atom

#endif // RNAMAKE_SRC_STRUCTURE_ALL_ATOM_ALIGNER_HPP_
