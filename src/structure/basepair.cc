//
// Created by Joseph Yesselman on 12/15/17.
//

#include <structure/basepair.h>

namespace structure {

primitives::BasepairType generate_bp_type(Residue const &res1,
                                          Residue const &res2,
                                          util::X3dnaBPType x3dna_bp_type) {

  if (x3dna_bp_type != util::X3dnaBPType::cWUW) {
    return primitives::BasepairType::NC;
  }

  auto bp_str = String();
  bp_str += res1.get_name();
  bp_str += res2.get_name();

  auto wc_names = Strings{"GC", "CG", "AU", "UA"};
  if (std::find(wc_names.begin(), wc_names.end(), bp_str) != wc_names.end()) {
    return primitives::BasepairType::WC;
  } else if (bp_str == "GU" || bp_str == "UG") {
    return primitives::BasepairType::GU;
  } else {
    return primitives::BasepairType::NC;
  };
}

} // namespace structure