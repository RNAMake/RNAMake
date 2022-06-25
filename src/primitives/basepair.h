//
// Created by Joseph Yesselman on 1/29/17.
//

#ifndef RNAMAKE_PRIMITIVES_BASEPAIR_H
#define RNAMAKE_PRIMITIVES_BASEPAIR_H

#include <util/uuid.h>
//#include <util/x3dna.h>
#include <primitives/residue.h>

/*
 * Exception for basepair
 */
class BasepairException : public std::runtime_error {
public:
  /**
   * Standard constructor for BasepairException
   * @param   message   Error message for basepair
   */
  explicit BasepairException(String const &message)
      : std::runtime_error(message) {}
};

namespace primitives {

// WC: watson crick basepair
// GU: gu basepair
// NC : mismatched basepair (non-conical)
enum class BasepairType { WC, GU, NC };

class Basepair {
public:
  [[nodiscard]] virtual const util::Uuid &
  get_partner(util::Uuid const &) const = 0;

  [[nodiscard]] virtual inline const BasepairType &get_bp_type() const = 0;

  [[nodiscard]] virtual inline util::Uuid const &get_uuid() const = 0;

  [[nodiscard]] virtual inline const String &get_name() const = 0;

  [[nodiscard]] virtual inline util::Uuid const &get_res1_uuid() const = 0;

  [[nodiscard]] virtual inline util::Uuid const &get_res2_uuid() const = 0;
};

template <typename Restype>
String generate_bp_name(Restype const &res1, Restype const &res2) {

  auto res1_name = String("");
  auto res2_name = String("");

  if (res1.get_i_code() == ' ') {
    res1_name = res1.get_chain_id() + std::to_string(res1.get_num());
  } else {
    res1_name = res1.get_chain_id() + std::to_string(res1.get_num()) +
                res1.get_i_code();
  }

  if (res2.get_i_code() == ' ') {
    res2_name = res2.get_chain_id() + std::to_string(res2.get_num());
  } else {
    res2_name = res2.get_chain_id() + std::to_string(res2.get_num()) +
                res2.get_i_code();
  }

  if (res1.get_chain_id() < res2.get_chain_id()) {
    return res1_name + "-" + res2_name;
  }
  if (res2.get_chain_id() < res1.get_chain_id()) {
    return res2_name + "-" + res1_name;
  }

  if (res1.get_num() < res2.get_num()) {
    return res1_name + "-" + res2_name;
  } else {
    return res2_name + "-" + res1_name;
  }
}

} // namespace primitives

#endif // TEST_BASEPAIR_H