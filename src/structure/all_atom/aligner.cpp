//
// Created by Joe Yesselman on 7/10/22.
//

#include <structure/all_atom/aligner.hpp>
// This is how you do alignment between two references
// after you get the data from db. Do it before adding to segment graph
namespace structure::all_atom {

void Aligner::align(const Segment &ref, Segment &seg, Index end_index) {

  _rot = math::rotation_between_frames(ref.get_end_ref_frame(end_index),
                                       seg.get_end_ref_frame(0));
  seg.rotate(_rot);
  seg.move(ref.get_end_center(end_index) - seg.get_end_center(0));

}

} // namespace structure::all_atom
