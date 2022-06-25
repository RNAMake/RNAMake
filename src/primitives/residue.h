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

public: // getters
  [[nodiscard]] virtual inline const util::Uuid &get_uuid() const = 0;
};

} // namespace primitives
#endif // PRIMITIVES_RESIDUE_H