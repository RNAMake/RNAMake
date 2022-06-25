//
// Created by Joseph Yesselman on 1/26/17.
//

#ifndef PRIMITIVES_RESIDUE_H
#define PRIMITIVES_RESIDUE_H

#include <base/types.hpp>
#include <util/uuid.h>

/*
 * Exception for residues
 */
class ResidueException : public std::runtime_error {
public:
  /**
   * Standard constructor for ResidueException
   * @param   message   Error message for residue
   */
  explicit ResidueException(const String &message)
      : std::runtime_error(message) {}
};

namespace primitives {

class Residue {
public:
  /**
   * @brief constructor for Residue class
   * @param   name        residue name (A, G, C, T)
   * @param   num         residue num
   * @param   chain_id    what chain does this residue belong to ("A", "B")
   * @param   uuid        residue unique indentifier
   */
  inline Residue(char name, int num, char chain_id, char i_code,
                 const util::Uuid &uuid)
      : _name(name), _num(num), _chain_id(chain_id), _i_code(i_code),
        _uuid(uuid) {}

  /**
   * copy construtor for residue class
   * @param  r   residue to be copied from
   */
  Residue(const Residue &r) = default;

  virtual ~Residue() = default;

protected:
  // let derived class setup members
  Residue() {}

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
  [[nodiscard]] virtual inline const util::Uuid &get_uuid() const = 0;
};

} // namespace primitives
#endif // PRIMITIVES_RESIDUE_H