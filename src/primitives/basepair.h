//
// Created by Joseph Yesselman on 1/29/17.
//

#ifndef RNAMAKE_PRIMITIVES_BASEPAIR_H
#define RNAMAKE_PRIMITIVES_BASEPAIR_H

#include <base/simple_string.h>
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
  BasepairException(String const &message) : std::runtime_error(message) {}
};

namespace primitives {

// WC: watson crick basepair
// GU: gu basepair
// NC : mismatched basepair (non-conical)
enum BasepairType { WC, GU, NC };

class Basepair {
public:
  inline Basepair(util::Uuid const &res1_uuid, util::Uuid const &res2_uuid,
                  util::Uuid const &uuid, BasepairType const &bp_type,
                  base::SimpleStringCOP const &name)
      : _res1_uuid(res1_uuid), _res2_uuid(res2_uuid), _uuid(uuid),
        _bp_type(bp_type), _name(name) {}

  inline Basepair(Basepair const &bp)
      : _res1_uuid(bp._res1_uuid), _res2_uuid(bp._res2_uuid), _uuid(bp._uuid),
        _bp_type(bp._bp_type), _name(bp._name) {}

  inline Basepair(util::Uuid const &res1_uuid, util::Uuid const &res2_uuid,
                  util::Uuid const &uuid)
      : _res1_uuid(res1_uuid), _res2_uuid(res2_uuid), _uuid(uuid),
        _bp_type(BasepairType::NC), _name(base::SimpleStringOP(nullptr)) {}

  virtual ~Basepair() {}

protected:
  Basepair() {}

public:
  /**
   * equal operator checks whether the unique indentifier is the same
   * @param other another basepair to check if its the same
   */

  inline bool operator==(Basepair const &other) const {
    return _uuid == other._uuid;
  }

  inline bool operator!=(Basepair const &other) const {
    return _uuid != other._uuid;
  }

public:
  util::Uuid const &get_partner(util::Uuid const &) const;

  inline BasepairType const &get_bp_type() const { return _bp_type; }

  inline util::Uuid const &get_uuid() const { return _uuid; }

  inline base::SimpleStringCOP get_name() const { return _name; }

  String get_name_str() const { return _name->get_str(); }

  inline util::Uuid const &get_res1_uuid() const { return _res1_uuid; }

  inline util::Uuid const &get_res2_uuid() const { return _res2_uuid; }

protected:
  util::Uuid _uuid;
  util::Uuid _res1_uuid, _res2_uuid;
  BasepairType _bp_type;
  base::SimpleStringCOP _name;
};

typedef Basepair PrimitiveBasepair;
typedef std::vector<PrimitiveBasepair> PrimitiveBasepairs;

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