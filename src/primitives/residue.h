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
  ResidueException(const String &message) : std::runtime_error(message) {}
};

namespace primitives {

class Residue {
public:
  /**
   * constructor for Residue class
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
  [[nodiscard]] inline util::Uuid const &get_uuid() const { return _uuid; }

  [[nodiscard]] virtual String get_str() const { return {""}; }

protected:
  char _name;
  int _num;
  char _chain_id;
  char _i_code;
  util::Uuid _uuid;
};

typedef Residue PrimitiveResidue;
typedef std::vector<Residue> PrimitiveResidues;

} // namespace primitives
#endif // PRIMITIVES_RESIDUE_H