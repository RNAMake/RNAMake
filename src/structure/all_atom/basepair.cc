//
// Created by Joseph Yesselman on 12/15/17.
//

#include <structure/all_atom/basepair.h>
#include <structure/all_atom/residue.h>

namespace structure::all_atom {

structure::base::BasepairType
generate_bp_type(Residue const &res1, Residue const &res2,
                 util::x3dna::X3dnaBPType x3dna_bp_type) {

  if (x3dna_bp_type != util::x3dna::X3dnaBPType::cWUW) {
    return structure::base::BasepairType::NC;
  }

  auto bp_str = String();
  bp_str += res1.get_name();
  bp_str += res2.get_name();

  auto wc_names = Strings{"GC", "CG", "AU", "UA"};
  if (std::find(wc_names.begin(), wc_names.end(), bp_str) != wc_names.end()) {
    return structure::base::BasepairType::WC;
  } else if (bp_str == "GU" || bp_str == "UG") {
    return structure::base::BasepairType::GU;
  } else {
    return structure::base::BasepairType::NC;
  };
}

} // namespace structure::all_atom