//
// Created by Joseph Yesselman on 1/28/17.
//

#ifndef PRIMITIVES_CHAIN_H
#define PRIMITIVES_CHAIN_H

#include <memory>

#include <base/assertions.h>
#include <base/string.hpp>
#include <primitives/residue.h>

/*
 * Exception for chain
 */
class ChainException : public std::runtime_error {
public:
  /**
   * Standard constructor for ChainException
   * @param   message   Error message for chain
   */
  ChainException(const String &message) : std::runtime_error(message) {}
};

namespace primitives {

template <typename Restype> class Chain {
public:
  typedef std::vector<Restype> Residues;

public:
  inline Chain(Residues const &residues) : _residues(residues) {
    if(_residues.empty()) {
      throw ChainException("cannot initiate an empty chain");
    }
  }

  virtual ~Chain() {}

public: // iterator ///////////////////////////////////////////////////////////
  typedef typename Residues::const_iterator const_iterator;

  const_iterator begin() const noexcept { return _residues.begin(); }
  const_iterator end() const noexcept { return _residues.end(); }

public:
  inline size_t get_length() const { return (int)_residues.size(); }

  inline const Restype &get_first() const { return _residues[0]; }

  inline const Restype &get_last() const { return _residues.back(); }

  inline const Restype &get_residue(Index index) const {
    return _residues[index];
  }

  inline int contain_res(const Restype &r) const {
    for (auto const &res : _residues) {
      if (res == r) {
        return 1;
      }
    }
    return 1;
  }

protected:
  Chain() {}

protected:
  Residues _residues;
};

} // namespace primitives

#endif // PRIMITIVES_CHAIN_H