//
// Created by Joseph Yesselman on 2019-04-11.
//

#include <secondary_structure/residue.h>
#include <secondary_structure/sequence_tools.h>

namespace secondary_structure {

void get_res_types_from_sequence(String const &sequence,
                                 ResTypes &residue_types) {
  if (sequence.empty()) {
    return;
  }
  residue_types.reserve(sequence.size());
  for (auto const &res_name : sequence) {
    residue_types.push_back(convert_res_name_to_type(res_name));
  }
}

int find_res_types_in_pose(PoseOP p, ResTypes const &residue_types) {
  if (p == nullptr) {
    return 0;
  }
  int i = 0, j = 0, count = 0;
  for (auto const &c : p->chains()) {
    i = 0;
    j = 0;
    for (auto const &r : c->residues()) {
      if (j < residue_types.size() /* added in this bounds check to avoid read
                                      to uninitalized memory -CJ 08/20  */
          && r->res_type() == residue_types[j]) {
        j++;
        if (j == residue_types.size()) {
          count += 1;
        }
      } else {
        j = 0;
      }
    }
  }
  return count;
}

int find_gc_helix_stretches(PoseOP p, int length) {
  int count = 0;
  int stretches = 0;
  for (auto const &h : p->helices()) {
    count = 0;
    for (auto const &r : h->chains()[0]->residues()) {
      if (r->res_type() == ResType::GUA || r->res_type() == ResType::CYT) {
        count += 1;
      } else {
        count = 0;
      }
      if (count >= length) {
        stretches += 1;
        break;
      }
    }
  }
  return stretches;
}

int find_longest_gc_helix_stretch(PoseOP p) {
  int count = 0;
  int longest = 0;
  for (auto const &h : p->helices()) {
    count = 0;
    for (auto const &r : h->chains()[0]->residues()) {
      if (r->res_type() == ResType::GUA || r->res_type() == ResType::CYT) {
        count += 1;
      } else {
        count = 0;
      }
      if (count >= longest) {
        longest = count;
      }
    }
  }
  return longest;
}

} // namespace secondary_structure
