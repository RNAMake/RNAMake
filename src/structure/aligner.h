//
// Created by Joseph Yesselman on 12/26/17.
//

#ifndef RNAMAKE_NEW_ALL_ATOM_ALIGNER_H
#define RNAMAKE_NEW_ALL_ATOM_ALIGNER_H

#include <base/log.h>
#include "primitives/aligner.h"
#include "structure/segment.h"
#include "structure/basepair.h"

namespace structure {

  class Aligner : public primitives::Aligner<Segment, Basepair>  {
  public:
      typedef primitives::Aligner<Segment, Basepair> BaseClass;
  public:
      Aligner():
              BaseClass() {
          dummy_ = math::Point();
      }

      ~Aligner() {}

  public:
      virtual
      void
      align(
              Basepair const & ref_bp,
              Segment & seg) const {

          r_ = dot(ref_bp.get_ref_frame().get_transposed(), seg.get_aligned_end().get_ref_frame());
          r_.unitarize();
          transpose(r_, r_t);
          t_ = -ref_bp.get_center();
          seg.transform(r_t, t_, dummy_);
          seg.move(ref_bp.get_center() - seg.get_aligned_end().get_center());

          sugar_dist_1_ = ref_bp.get_res1_c1_prime_coord().distance(seg.get_aligned_end().get_res1_c1_prime_coord());
          sugar_dist_2_ = ref_bp.get_res2_c1_prime_coord().distance(seg.get_aligned_end().get_res1_c1_prime_coord());
          if(sugar_dist_1_ > 5 && sugar_dist_2_ > 5) {
              LOGW << "difference in sugar c1' coords between reference and aligned is greater than 5. " <<
                   "This could lead to alignment issues!";
              return;
          }

          if(sugar_dist_1_ < sugar_dist_2_) {
              sugar_diff_1_ = ref_bp.get_res1_c1_prime_coord() - seg.get_aligned_end().get_res1_c1_prime_coord();
              sugar_diff_2_ = ref_bp.get_res2_c1_prime_coord() - seg.get_aligned_end().get_res2_c1_prime_coord();
          }
          else {
              sugar_diff_1_ = ref_bp.get_res1_c1_prime_coord() - seg.get_aligned_end().get_res2_c1_prime_coord();
              sugar_diff_2_ = ref_bp.get_res2_c1_prime_coord() - seg.get_aligned_end().get_res1_c1_prime_coord();
          }

          avg_sugar_diff_ = (sugar_diff_1_ + sugar_diff_2_) / 2;
          seg.move(avg_sugar_diff_);

      }

      virtual
      SegmentTypeOP
      get_aligned(
              Basepair const & ref_bp,
              Segment const & seg) const {
          auto seg_copy = std::make_shared<Segment>(seg);
          align(ref_bp, *seg_copy);
          return seg_copy;
      }

  private:
      mutable math::Matrix r_, r_t;
      mutable math::Point t_, dummy_, sugar_diff_1_, sugar_diff_2_, avg_sugar_diff_;
      mutable double sugar_dist_1_, sugar_dist_2_;

  };

}


#endif //RNAMAKE_NEW_ALIGNER_H