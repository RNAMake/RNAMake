//
// Created by Joe Yesselman on 6/22/22.
//

#include "../common.hpp"

#include <primitives/basepair.h>
#include <primitives/chain.h>
#include <primitives/residue.h>

// try trival implementations ///////////////////////////////////////////////

class Residue : public primitives::Residue {
public:
  inline Residue(char name, int num, char chain_id, char i_code,
                 const util::Uuid &uuid)
      : _name(name), _num(num), _chain_id(chain_id), _i_code(i_code),
        _uuid(uuid) {}

  Residue(const Residue &r) = default;

  virtual ~Residue() = default;

public:
  /**
   * equal operator checks whether the unique indentifier is the same
   * @param   r   another residue to check if its the same
   */
  inline bool operator==(const Residue &other) const {
    return _uuid == other._uuid;
  }

  inline bool operator!=(const Residue &other) const {
    return _uuid != other._uuid;
  }

public: // getters
  /**
   * getter the chain_id, i.e. "A", "B", the id of the chain this residue
   * belongs to
   */
  [[nodiscard]] inline char get_chain_id() const { return _chain_id; }

  /**
   * getter for the name of the residue, i.e. "A", "G" etc
   */
  [[nodiscard]] inline char get_name() const { return _name; }

  /**
   * getter for the residue num
   */
  [[nodiscard]] inline int get_num() const { return _num; }

  /**
   * getter for the residue insertion code
   */
  [[nodiscard]] inline char get_i_code() const { return _i_code; }

  /**
   * getter for residue unique indentifier
   */
  [[nodiscard]] inline util::Uuid const &get_uuid() const override {
    return _uuid;
  }

private:
  char _name;
  int _num;
  char _chain_id;
  char _i_code;
  util::Uuid _uuid;
};

class Basepair : public primitives::Basepair {
public:
  inline Basepair(const util::Uuid &res1_uuid, const util::Uuid &res2_uuid,
                  const util::Uuid &uuid,
                  const primitives::BasepairType &bp_type, String &name)
      : _res1_uuid(res1_uuid), _res2_uuid(res2_uuid), _uuid(uuid),
        _bp_type(bp_type), _name(std::move(name)) {}

  inline Basepair(Basepair const &bp) = default;

  virtual ~Basepair() = default;

public:
  inline bool operator==(Basepair const &other) const {
    return _uuid == other._uuid;
  }

  inline bool operator!=(Basepair const &other) const {
    return _uuid != other._uuid;
  }

public:
  [[nodiscard]] util::Uuid const &
  get_partner(util::Uuid const &uuid) const override {
    if (uuid == _res1_uuid) {
      return _res2_uuid;
    } else {
      return _res1_uuid;
    }
  }

  [[nodiscard]] inline primitives::BasepairType const &
  get_bp_type() const override {
    return _bp_type;
  }

  [[nodiscard]] inline util::Uuid const &get_uuid() const override {
    return _uuid;
  }

  [[nodiscard]] inline const String &get_name() const override { return _name; }

  [[nodiscard]] inline util::Uuid const &get_res1_uuid() const override {
    return _res1_uuid;
  }

  [[nodiscard]] inline util::Uuid const &get_res2_uuid() const override {
    return _res2_uuid;
  }

private:
  util::Uuid _uuid;
  util::Uuid _res1_uuid;
  util::Uuid _res2_uuid;
  primitives::BasepairType _bp_type;
  String _name;
};


TEST_CASE("brief test of primitive functionility") {
  SUBCASE("test residue") {
    Residue r1 = {'A', 1, 'A', ' ', util::generate_uuid()};
    CHECK(r1.get_name() == 'A');
    CHECK(r1.get_num() == 1);
    auto uuid = r1.get_uuid();
    CHECK(r1.get_uuid() == uuid);
    Residue r2 = {'A', 1, 'A', ' ', uuid};
    CHECK(r1 == r2);
  }
  SUBCASE("test basepairs") {
    SUBCASE("test default construction") {
      String name = "A1-A2";
      Basepair bp1 = {util::generate_uuid(), util::generate_uuid(),
                      util::generate_uuid(), primitives::BasepairType::WC, name};
      CHECK(bp1.get_name() == "A1-A2");
    }
    SUBCASE("test name construction") {
      Residue r1 = {'A', 1, 'A', ' ', util::generate_uuid()};
      Residue r2 = {'A', 2, 'A', ' ', util::generate_uuid()};
      String bp_name = primitives::generate_bp_name<Residue>(r1, r2);
      CHECK(bp_name == "A1-A2");
    }
  } /*
  SUBCASE("test chain") {
    Residue r1 = {'A', 1, 'A', ' ', util::generate_uuid()};
    Residue r2 = {'A', 2, 'A', ' ', util::generate_uuid()};
    Chain<Residue> c({r1, r2});
  }*/
}