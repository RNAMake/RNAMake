//
// Created by Joseph Yesselman on 12/26/17.
//

#ifndef RNAMAKE_NEW_ALIGNER_H
#define RNAMAKE_NEW_ALIGNER_H

#include <memory>

namespace primitives {

template <typename SegmentType, typename BasepairType> class Aligner {
public:
  typedef std::shared_ptr<SegmentType> SegmentTypeOP;

public:
  Aligner() {}

  ~Aligner() {}

public:
  virtual void align(BasepairType const &ref_bp, SegmentType &seg) const = 0;

  virtual SegmentTypeOP get_aligned(BasepairType const &ref_bp,
                                    SegmentType const &seg) const = 0;
};

} // namespace primitives

#endif // RNAMAKE_NEW_ALIGNER_H